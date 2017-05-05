# -*- coding: utf-8 -*-
# arg1 : toysequence.fasta


########################## trouver un moyen pour mettre et changer le threshold !


################################ IMPORT ################################################################################################

#from sequana import sequence
import subprocess
import pandas as pd
import sys
import os

################################ PARAMETERS ##############################################################################################

file_sequence =  str(sys.argv[1])


################################ INPUT DATA ################################################################################################

#seq = sequence.Sequence(file_sequence)

## example

# df = pd.DataFrame()
# list_df = []
# task = subprocess.Popen(["cat", file_sequence], stdout=subprocess.PIPE)
# i = 0
# for line in task.stdout:
#     i +=1
#     list_df.append(line)
#     #df=pd.concat([df, line])
# task.wait()



# # get first line of fasta (needed in input shustring) and replace spaces
# first_line = subprocess.check_output(["head","-n","1",file_sequence])
# first_line = first_line.decode('utf8').replace("\n","").replace(" ","_")

# # read fasta
# task_read           = subprocess.Popen(["cat", file_sequence], stdout=subprocess.PIPE)
# # replace spaces of fasta
# task_replace_spaces = subprocess.Popen(["sed","s/ /_/g"], stdin=task_read.stdout, stdout=subprocess.PIPE)
# # shustring command
# task_shus           = subprocess.Popen(['shustring','-r','-q','-l',first_line], stdin=task_replace_spaces.stdout, stdout=subprocess.PIPE)
# # read stdout line by line and append to list
# list_df = []
# for line in task_shus.stdout:
#     list_df.append(line.decode('utf8').replace("\n",'').split("\t"))
#     #df=pd.concat([df, line])
# task_shus.wait()
# df=pd.DataFrame(list_df[2:len(list_df)])

# print(list_df)
# print(df)

# # the first line of shustring output contains genome length, longest shustring
# len_genome        = int(list_df[0][1])
# longest_shustring = int(list_df[0][3].split("<=")[2])






class Repeats(object):
    """Class for finding repeats in DNA or RNA linear sequences.

    .. note:: You may use a Fasta file as input


    """
    def __init__(self, filename_fasta):
        """.. rubric:: Constructor

        Input must be a fasta file with valid DNA or RNA characters

        :param str filename_fasta: String of a Fasta File, only the first
         sequence is used.
        :param int threshold: Minimal length of repeat to output
        """
        self.threshold                  = None
        self._df_shustring              = None
        self._header                    = None
        self._length                    = None
        self._longest_shustring         = None
        self._begin_end_repeat_position = None
        self._filename_fasta            = filename_fasta
        self._previous_thr              = None

    def _get_header(self):
        """get first line of fasta (needed in input shustring) 
        and replace spaces by underscores
        """
        if self._header is None:
            first_line = subprocess.check_output(["head","-n","1",file_sequence])
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


    def _find_begin_end_repeats(self):
        """Returns position of repeats longer thant threshold
        as an ordered dataframe
        """
        # if a threshold is already set, save the previous one
        if self.threshold is not None:
            previous_thr = self.threshold
        else:
            previous_thr = None
        
        self.threshold = threshold

        if self.df_shustring is None:
            self._get_shustrings_length()

        # if there is n result, or the threshold has changed
        if (self._begin_end_repeat_position is None) | (self.threshold != previous_thr):
            nb_row = self.df_shustring.shape[0]
            i = 0
            step_repeat_seq = []
            be_repeats = []
            e = 0

            # use list because much faster
            list_len_shus = list(self.df_shustring.loc[:,"shustring_length"])

            while(i < nb_row):
                # begining of repeat
                if (list_len_shus[i] > self.threshold):
                    b = i
                    # add 0 for non repeat area
                    step_repeat_seq = step_repeat_seq + [0 for j in range(e,b,1)]
                    # compute new end of repeat
                    len_repeat = list_len_shus[i]
                    e = b + len_repeat
                    # add len of repeat for all repeat area
                    step_repeat_seq = step_repeat_seq + [len_repeat for j in range(b,e,1)]
                    # save (b,e)
                    be_repeats.append((b,e))
                    # update i
                    i = e-1
                i +=1
            # add 0 at the end if needed
            step_repeat_seq = step_repeat_seq + [0 for j in range(e,i,1)]

            # if there are repeats, merge repeats that are fragmented
            if len(be_repeats) > 0:
                prev_tup = be_repeats[0]
                b = prev_tup[0]
                begin_end_repeat_position = []
                for tup in be_repeats[1:len(be_repeats)]:
                    if tup[0] == prev_tup[1]:
                        # concat
                        e = tup[1]
                        prev_tup = tup
                    else:
                        # real end of repeat : append result and update b, e
                        e = prev_tup[1]
                        begin_end_repeat_position.append((b,e))
                        prev_tup = tup
                        b = prev_tup[0]
            else:
                begin_end_repeat_position = []
            self._begin_end_repeat_position = begin_end_repeat_position

        return self._begin_end_repeat_position
    begin_end_repeat_position = property(_find_begin_end_repeats)








