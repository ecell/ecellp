import sys
sys.path.append('../')

from generators import Genbank
from generators import Promoter
from generators import Terminator

def generator_test():
    conf = '../../../conf.ini'
    gb = Genbank(conf)
    gb.generate_genbank_file()
    print "Finished gbk file"

    tm = Terminator(conf)
    tm.generate_terminator_file()
    print "Finished terminator file"
    
    pr = Promoter(conf)
    pr.generate_promoter_file()
    print "Finisehd promoter file"

if __name__ == '__main__':
    generator_test()



