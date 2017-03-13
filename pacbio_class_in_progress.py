from sequana.lazy import pylab
from sequana.lazy import numpy as np
from sequana.lazy import pandas as pd
import collections #lazy ?
import pysam
from biokit.viz import hist2d



class PacBioBAM(object):
    """
    """
    def __init__(self, filename):
        self.filename = filename
        self.data = pysam.AlignmentFile(filename, check_sq=False)
        self._N = None
        self._df = None

    def __len__(self):
        if self._N is None:
            df = self._get_df()
        return self._N

    def reset(self):
        self.data.close()
        self.data = pysam.AlignmentFile(self.filename, check_sq=False)

    def stride(self, output_filename, stride=10):
        self.reset()
        with pysam.AlignmentFile(output_filename,  "wb", template=self.data) as fh:

            for i, read in enumerate(self.data):
                if i % stride == 0: 
                    fh.write(read)
                    
    def _get_df(self):
        if self._df is None:
            self.reset()
            N = 0
            
            all_results = []
            for read in self.data:
                res = []
                # count reads
                N += 1
                if (N % 10000) == 0:
                    print("Read %d sequences" %N)
                #res[0] = read length
                res.append(read.query_length)
                # res[1] = GC content
                c = collections.Counter(read.query_sequence)
                res.append( (c['g'] + c['G'] + c['c'] + c['C'])/float(sum(c.values())) )
                # res[2] = snr A
                # res[3] = snr C
                # res[4] = snr G
                # res[5] = snr T
                snr = list([x for x in read.tags if x[0]=='sn'][0][1])
                res = res + snr
                #res[6] = ZMW name
                res.append(read.qname.split('/')[1])
                
                # aggregate results
                all_results.append(res)

            self._df = pd.DataFrame(all_results, columns=['read_length','GC_content','snr_A','snr_C','snr_G','snr_T','ZMW'])
            self._N = N
            self.reset()     
        return self._df
