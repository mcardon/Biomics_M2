# -*- coding: utf-8 -*-
# arg1 : toysequence.fasta


########################## trouver un moyen pour mettre et changer le threshold !


################################ IMPORT ################################################################################################

#from sequana import sequence
import subprocess
import pandas as pd
import sys
import pylab

################################ PARAMETERS ##############################################################################################

file_sequence =  str(sys.argv[1])


################################ INPUT DATA ################################################################################################


class Repeats(object):
    """Class for finding repeats in DNA or RNA linear sequences.

    .. note:: You may use a Fasta file as input


    """
    def __init__(self, filename_fasta,merge=False):
        """.. rubric:: Constructor

        Input must be a fasta file with valid DNA or RNA characters

        :param str filename_fasta: String of a Fasta File, only the first
         sequence is used.
        :param int threshold: Minimal length of repeat to output
        """
        self._threshold                       = None
        self._df_shustring                    = None
        self._header                          = None
        self._length                          = None
        self._longest_shustring               = None
        self._begin_end_repeat_position       = None
        self._begin_end_repeat_position_merge = None
        self._filename_fasta                  = filename_fasta
        self._previous_thr                    = None
        self._list_len_repeats                = None
        if not isinstance(merge, bool):
            raise TypeError("do_merge must be boolean")
        self._do_merge                        = merge

    def _get_header(self):
        """get first line of fasta (needed in input shustring) 
        and replace spaces by underscores
        """
        if self._header is None:
            first_line = subprocess.check_output(["head","-n","1",self._filename_fasta])
            first_line = first_line.decode('utf8')
            first_line = first_line.replace("\n","").replace(" ","_")
            self._header = first_line
        return self._header
    header = property(_get_header)

    def _get_shustrings_length(self):
        """Return dataframe with shortest unique substring length at each position
        shortest unique substrings are unique in the sequence and its complement
        Uses shustring tool"""
        if self._df_shustring is None:
            # read fasta
            task_read           = subprocess.Popen(["cat", self._filename_fasta], 
                stdout=subprocess.PIPE)
            # replace spaces of fasta
            task_replace_spaces = subprocess.Popen(["sed","s/ /_/g"], 
                stdin=task_read.stdout, stdout=subprocess.PIPE)
            # shustring command
            task_shus           = subprocess.Popen(['shustring','-r','-q','-l',self.header], 
                stdin=task_replace_spaces.stdout, stdout=subprocess.PIPE)

            # read stdout line by line and append to list
            list_df = []
            for line in task_shus.stdout:
                list_df.append(line.decode('utf8').replace("\n",'').split("\t"))
                #df=pd.concat([df, line])
            task_shus.wait()

            # convert to dataframe
            df = pd.DataFrame(list_df[2:len(list_df)])
            self._df_shustring = df.astype(int)
            self._df_shustring.columns = ["position","shustring_length"]

            # get input sequence length and longest shustring in the first line
            self._length = int(list_df[0][1])
            self._longest_shustring = int(list_df[0][3].split("<=")[2])

        return self._df_shustring
    df_shustring = property(_get_shustrings_length)

    def _get_genome_length(self):
        if self._df_shustring is None:
            self._get_shustrings_length()
        return self._length
    length = property(_get_genome_length)

    def _get_longest_shustring(self):
        if self._df_shustring is None:
            self._get_shustrings_length()
        return self._longest_shustring
    longest_shustring = property(_get_longest_shustring)


    def _find_begin_end_repeats(self,force=False):
        """Returns position of repeats longer than threshold
        as an ordered list
        """
        if self.df_shustring is None:
            self._get_shustrings_length()

        if self._threshold is None:
            #print("No threshold : please set minimul length of repeats to output")
            raise ValueError("threshold : please set threshold (minimum length of repeats to output)")

        # if there is no result yet, or the threshold has changed
        if (self._begin_end_repeat_position is None) | (self.threshold != self._previous_thr) | force:
            nb_row = self.df_shustring.shape[0]
            i = 0
            step_repeat_seq = []
            be_repeats = []
            e = 0

            # use list because faster
            list_len_shus = list(self.df_shustring.loc[:,"shustring_length"])

            while(i < nb_row):
                # begining of repeat
                if (list_len_shus[i] > self.threshold):
                    b = i
                    # compute new end of repeat
                    len_repeat = list_len_shus[i]
                    e = b + len_repeat
                    # save (b,e)
                    be_repeats.append((b,e))
                    # update i
                    i = e-1
                i +=1

            self._begin_end_repeat_position = be_repeats

        self._get_merge_repeats()

            

        #return self._begin_end_repeat_position
    def _get_be_repeats(self):
        self._find_begin_end_repeats()
        return self._begin_end_repeat_position

    begin_end_repeat_position = property(_get_be_repeats)

    def _set_threshold(self,value):
        if value < 0:
            raise ValueError("Threshold must be positive integer")
        if value != self._threshold:
            self._previous_thr = self._threshold
        self._threshold = value
        self._find_begin_end_repeats()
        self._list_len_repeats = [tup[1]-tup[0] for tup in self._begin_end_repeat_position]

    def _get_threshold(self):
        return self._threshold

    threshold = property(_get_threshold, _set_threshold)

    def _get_list_len_repeats(self):
        if self._list_len_repeats is None:
            raise UserWarning("Please set threshold (minimum length of repeats to output)")
        return self._list_len_repeats

    list_len_repeats = property(_get_list_len_repeats)

    def _get_merge_repeats(self):
        if self._do_merge:
            # if there are repeats, merge repeats that are fragmented
            if len(self._begin_end_repeat_position) > 0:
                prev_tup = self._begin_end_repeat_position[0]
                b = prev_tup[0]
                begin_end_repeat_position_merge = []
                for i in range(1,len(self._begin_end_repeat_position),1):
                    tup = self._begin_end_repeat_position[i]
                    if tup[0] == prev_tup[1]:
                        # concat
                        e = tup[1]
                        prev_tup = tup
                        if i == (len(self._begin_end_repeat_position) -1):
                            # last tup : append to result
                            begin_end_repeat_position_merge.append((b,e))
                    else:
                        # real end of repeat : append result and update b, e
                        e = prev_tup[1]
                        begin_end_repeat_position_merge.append((b,e))
                        prev_tup = tup
                        b = prev_tup[0]

            self._begin_end_repeat_position = begin_end_repeat_position_merge

    def _get_do_merge(self):
        return self._do_merge

    def _set_do_merge(self, do_merge):
        if not isinstance(do_merge, bool):
            raise TypeError("do_merge must be boolean")
        # if different
        if do_merge != self._do_merge:
            self._do_merge = do_merge
            if self._do_merge:
                # did not merge before, merge now
                if self._begin_end_repeat_position is None:
                    self._find_begin_end_repeats()
            else:
                # data is already merged : need to compute again to un-merge
                self._find_begin_end_repeats(force=True)

    do_merge = property(_get_do_merge,_set_do_merge)

    def hist_length_repeats(self, bins=50, alpha=0.5, hold=False, fontsize=12,
                        grid=True, label="",xlabel="Repeat length", ylabel="#"):
        # check that user has set a threshold
        if self._list_len_repeats is None:
            self._get_list_len_repeats()

        if hold is False:
            pylab.clf()
        pylab.hist(self._list_len_repeats, alpha=alpha, label=label, bins=bins)
        pylab.legend()
        pylab.xlabel(xlabel, fontsize=fontsize)
        pylab.ylabel(ylabel, fontsize=fontsize)
        if grid is True:
            pylab.grid(True)








