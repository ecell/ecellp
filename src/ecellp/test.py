__program__ = 'test'
__version__ = '0.1'
__author__ = 'Kazunari Kaizu <kaizu@riken.jp>'
__copyright__ = ''
__license__ = ''

import core
import domain
import sequence_utils


def test1(filename):
    algn = core.Alignment(
        sequence_utils.read_sequence(filename), name = "pTAK117", type = "PLASMID")
    print sequence_utils.format_sequence(algn.sequence())
    print algn.is_cyclic()
    print algn.size()
    print algn.annotations()
    print algn.sequence(1, 10, +1)
    print algn.sequence(10, 20, +1)
    print algn.sequence(6080, 10, +1)
    print algn.sequence(20, 10, -1)
    print algn.sequence(10, 6080, -1)

def test2(filename):
    algn = core.Alignment(
        sequence_utils.read_sequence(filename), name = "pTAK117", type = "PLASMID")

    print "Ptrc", algn.sequence(11, 98, +1)
    print "PL-s1con", algn.sequence(6067, 5496, -1)
    print "CI", algn.sequence(99, 714, +1)
    print "GFP", algn.sequence(914, 1630, +1)
    print "LacI", algn.sequence(5489, 4407, -1)
    # print "rbs E", algn.sequence(87, 101, +1)
    # print "rbs B", algn.sequence(894, 907, +1)
    print "rrn T1", algn.sequence(1846, 1889, +1)
    print "rrn T2", algn.sequence(2021, 2048, +1)
    print "rrn T1", algn.sequence(-1903, -1946, -1)
    print "rrn T2", algn.sequence(-2078, -2105, -1)

    print "TSS1", algn.sequence(60, 60, +1)
    print "TSS2", algn.sequence(-94, -94, -1)

    print "Olac", algn.sequence(60, 81, +1)
    print "OL3", algn.sequence(6061, 6045, -1)
    print "OL2", algn.sequence(6040, 6025, -1)
    print "OL1", algn.sequence(6016, 6001, -1)

    # region = algn.find("ACCACTGGCGGTTATA", -1)
    # print region, algn.sequence(*region)
    # print "", algn.sequence(24, 53, +1)
    # print "", algn.sequence(6027, 5999, -1)

def test3(filename):
    algn = core.Alignment(
        sequence_utils.read_sequence(filename), name = "pTAK117", type = "PLASMID")
    trk = core.Track(algn.size(), algn.is_cyclic())

    trk.place_domain(domain.Domain("A"), 1)
    trk.place_domain(domain.Domain("B"), 2)
    trk.place_domain(domain.Domain("A"), 3)
    print trk.list_domains()
    print trk.states()[: 10]
    print "removed", trk.remove_domain(1)
    print trk.list_domains()
    print trk.states()[: 10]
    print "removed", trk.remove_domain(2)
    print trk.list_domains()
    print trk.states()[: 10]
    trk.place_domain(domain.Domain("C"), 1)
    trk.place_domain(domain.Domain("C"), 2)
    print trk.list_domains()
    print trk.states()[: 10]
    print trk.query_domains_by_region(1, 2, +1)
    print trk.query_domains_by_region(1, 3, +1)
    trk.place_domain(domain.Domain("D"), 4, 8, +1)
    print trk.states()[: 10], trk.list_domains()

def test4(filename):
    algn = core.Alignment(
        core.read_sequence(filename), name = "pTAK117", type = "PLASMID")

    algn.place_domain(domain.Domain("Ptrc", type = "PROMOTER"), 11, 98, +1)
    algn.place_domain(domain.Domain("PL-s1con", type = "PROMOTER"), 6067, 5496, -1)
    algn.place_domain(domain.Domain("CI", type = "GENE"), 99, 714, +1)
    algn.place_domain(domain.Domain("GFP", type = "GENE"), 914, 1630, +1)
    algn.place_domain(domain.Domain("LacI", type = "GENE"), 5489, 4407, -1)
    algn.place_domain(domain.Domain("rrn T1", type = "TERMINATOR"), 1846, 1889, +1)
    algn.place_domain(domain.Domain("rrn T2", type = "TERMINATOR"), 2021, 2048, +1)
    algn.place_domain(domain.Domain("rrn T1", type = "TERMINATOR"), -1903, -1946, -1)
    algn.place_domain(domain.Domain("rrn T2", type = "TERMINATOR"), -2078, -2105, -1)
    algn.place_domain(domain.Domain("TSS", type = "TSS"), 60, 60, +1) # TSS1
    algn.place_domain(domain.Domain("TSS", type = "TSS"), -94, -94, -1) # TSS2
    algn.place_domain(domain.Domain("Olac", type = "OPERATOR"), 60, 81, +1)
    algn.place_domain(domain.Domain("OL3", type = "OPERATOR"), 6061, 6045, -1)
    algn.place_domain(domain.Domain("OL2", type = "OPERATOR"), 6040, 6025, -1)
    algn.place_domain(domain.Domain("OL1", type = "OPERATOR"), 6016, 6001, -1)
    algn.place_domain(domain.Domain("SFBS", type = "SFBS"), 24, 53, +1)
    algn.place_domain(domain.Domain("SFBS", type = "SFBS"), 6027, 5999, -1)
    print algn.num_domains()

    doms = algn.query_domains_by_region(-50, 50, +1)
    print doms
    print filter(lambda dom: dom.get_attribute("stride") > 0, doms)
    print filter(
        lambda dom: dom.get_attribute("type") == "TERMINATOR",
        algn.query_domains_by_region(-1904, +1847, +1))

    return algn

def test5(filename):
    pool = core.CompartmentSpace(1e-18)
    print pool.volume()
    pool.add_molecules(domain.Domain("A"), 10)
    pool.add_molecules(domain.Domain("B"), 30)
    print pool.list_domains()
    print pool.num_molecules(domain.Domain("A"))
    print pool.num_molecules(domain.Domain("B"))
    pool.remove_molecules(domain.Domain("B"), 15)
    pool.remove_molecules(domain.Domain("A"), 10)
    print pool.list_domains()
    print pool.num_molecules(domain.Domain("A"))
    print pool.num_molecules(domain.Domain("B"))

def test6():
    volume = 1e-18
    comp = core.CompartmentSpace(volume)
    comp.add_molecules(domain.Domain("RNAP", type = "MOLECULE"), 2000)
    comp.add_molecules(domain.Domain("SIGMA", type = "MOLECULE"), 1000)
    return comp

def test7(algn, comp):
    trk = algn.create_track("Binding proteins")
    trk.place_domain(domain.Domain("SIGMA", type = "MOLECULE"), 50)

    cdom1 = domain.ConditionalDomain("TSS")
    cdom1.add_condition(
        domain.Condition(-40, 0, includes = [domain.Domain("SIGMA")]))
    rr1 = core.SecondOrderInteraction(cdom1, domain.Domain("RNAP"), None, 1.0)

    cdom1 = domain.ConditionalDomain("SFBS")
    cdom1.add_condition(
        domain.Condition(
            0, +29, excludes = [domain.FilteringDomain(type = "MOLECULE")]))
    rr2 = core.SecondOrderInteraction(cdom1, domain.Domain("SIGMA"), None, 1.0)

    for rr in (rr1, rr2, ):
        r1, r2 = rr.reactants()
        print r1.find(algn), comp.num_molecules(r2)


if __name__ == '__main__':
    import sys


    test1(sys.argv[1])
    test2(sys.argv[1])
    test3(sys.argv[1])
    algn = test4(sys.argv[1])
    test5(sys.argv[1])
    comp = test6()
    test7(core.AlignmentSpace(algn), comp)
