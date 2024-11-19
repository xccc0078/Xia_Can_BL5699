import rna_draw as rd
import nupack
from nupack import *
from Bio.Seq import Seq

def draw_complex(result,complex):
    seqstr = complex.seq()
    dotpar = str(result[complex].mfe[0].structure)
    return rd.rna_draw(ss=dotpar,
                       seq=seqstr,
                render_type='strand',color_palette='pastel')



def saber_plish_probe(b1,b2,primer,probe_len=42,pos=0,toehold=False,toe_len=4,stem_len=25):
    """return a strand for a saber plish probe, based on the input bridges
    b1: immuno_saber bridge
    b2: saber_fish bridge
    primer: primer sequence
    probe_length: desired saber_plish probe_length
    pos: relative position from centre. 0 means equal, negative means more hybridisation on immunosaber,
    positive means more hybridisation on saber-fish
    """
    is_len = (probe_len // 2) - pos # ImmunoSaber part
    sf_len = probe_len - is_len # SaberFish part
    target = str(b1)[len(b1)-is_len:]+ str(b2)[:sf_len]
    probe = str(Seq(target).reverse_complement())
    assert len(probe) == probe_len, f'len(probe) = {len(probe)} for {pos} with sf_len={sf_len}'
    name = f"saber-plish({pos})"
    if toehold:
        assert((toe_len + stem_len) < probe_len)
        toe = probe[:toe_len]
        s1 = probe[toe_len:(toe_len+stem_len)]
        l1 = probe[(toe_len+stem_len):probe_len]
        s2 = str(Seq(s1).reverse_complement())
        probe = toe + s1 + l1 + s2
        name= f"{name}-th({toe_len},{stem_len})"
    return nupack.Strand(probe + "TT" + primer,name=name)


def report_complexes(IS_bridge,SF_bridge,primer,probe_len=42,pos=0,
                     sodium=.390,celsius=42,
                     toehold=False,toe_len=4,stem_len=25):
                         
    SP_probe = saber_plish_probe(IS_bridge,SF_bridge,primer,probe_len=probe_len,pos=pos,
                     toehold=toehold,toe_len=toe_len,stem_len=stem_len)

    ligation_bridge = Strand(str(IS_bridge) + str(SF_bridge),name='Ligation bridge')

    c3 = Complex([SF_bridge,SP_probe,IS_bridge])
    c3L = Complex([ligation_bridge,SP_probe]) 
    bc1 = Complex([SF_bridge,SP_probe])
    bc2 = Complex([IS_bridge,SP_probe])

    tube_model = Model(material='dna',sodium=sodium,celsius=celsius)
    
    saberplish_unligated = Tube(strands={SF_bridge: 1e-6,IS_bridge: 1e-6,SP_probe: 1e-6},
                                complexes=SetSpec(max_size=3,
                                                  include=[c3,bc1,bc2]),
                                name='unligated SABER-PLISH')
                         
    saberplish_ligated = Tube(strands={ligation_bridge: 1e-6,SP_probe: 1e-6},
                                complexes=SetSpec(max_size=2,
                                                  include=[c3L]),
                                name='ligated SABER-PLISH')

    saberplish_no_is = Tube(strands={SF_bridge: 1e-6,SP_probe: 1e-6},
                                complexes=SetSpec(max_size=2,
                                                  include=[bc1]),
                                name='neg control (no-immunoSaber)')
    saberplish_no_sf = Tube(strands={IS_bridge: 1e-6, SP_probe: 1e-6},
                                complexes=SetSpec(max_size=2,
                                                  include=[bc2]),
                                name='negative control (no saber-FISH)')

# run tube analysis job
    TA = tube_analysis(tubes=[saberplish_unligated,saberplish_ligated,saberplish_no_is,saberplish_no_sf],
                        model=tube_model,
                        compute=['pfunc','mfe'])
    return {'strand': str(SP_probe),
            'pos' : pos,
            'probe_len': probe_len,
            'celsius' : celsius,
            'sodium' : sodium,
            'toehold' : toehold,
            'toe_len' : toe_len,
            'stem_len' : stem_len,
            'c3':     TA[saberplish_unligated].complex_concentrations[c3],
            'c3L':     TA[saberplish_ligated].complex_concentrations[c3L],
            'bc1': TA[saberplish_unligated].complex_concentrations[bc1],
            'bc2': TA[saberplish_unligated].complex_concentrations[bc2],
            'bc1_nc': TA[saberplish_no_is].complex_concentrations[bc1],
            'bc2_nc': TA[saberplish_no_sf].complex_concentrations[bc2]}
                         

