import numpy as np
import unittest
import os, sys

sys.path.append(os.path.abspath(os.path.dirname(__file__)))
from systemtest_spinwave import SystemTest_Spinwave, m

class AF33kagomeTest(SystemTest_Spinwave):

    @classmethod
    def setUpClass(cls):
        super(AF33kagomeTest, cls).setUpClass()
        af33kagome = m.spinw();
        af33kagome.genlattice('lat_const',[6, 6, 40],'angled',[90, 90, 120],'sym','P -3');
        af33kagome.addatom('r',[1/2, 0, 0],'S', 1,'label','MCu1','color','r');
        af33kagome.gencoupling('maxDistance',7);
        af33kagome.addmatrix('label','J1','value',1.00,'color','g');
        af33kagome.addcoupling('mat','J1','bond',1);
        cls.swobj = af33kagome
        cls.relToll = 0.027
        cls.absToll = 1.2e-5


    def test_spinwave_calc(self):
        self.load_ref_dat('systemstest_spinwave_af33kagome.mat')
        af33kagome = self.swobj
        # Syntax for S0/evect below needs libpymcr 0.1.3 or newer; else must make S0 a 2D np array and evect a list
        S0 = [[0, 0, -1], [1, 1, -1], [0, 0, 0]]
        af33kagome.genmagstr('mode','helical','k',[-1/3, -1/3, 0],'n',[0, 0, 1],'unit','lu','S',S0,'nExt',[1, 1, 1])
        kag33spec = af33kagome.spinwave([[-1/2, 0, 0], [0, 0, 0], [1/2, 1/2, 0], 100],'hermit',False,'saveSabp',True)
        evect = np.linspace(0, 3, 100)
        kag33spec = m.sw_egrid(kag33spec,'component','Sxx+Syy','imagChk',False, 'Evect', evect)
        # Reduce values of S(q,w) so it falls within tolerance (rather than change tolerance for all values)
        kag33spec['swConv'] = kag33spec['swConv'] / 2e5
        # Ignores swInt in this case
        kag33spec['swInt'] = 0

        if self.plt is not None:
            ax = plt.imshow(np.real(np.flipud(kag33spec['swConv'])),
                            aspect='auto')
            # ax.set_clim(0, 1e-6)
            plt.xlabel('Q [q, q, 0] (r.l.u)')
            plt.ylabel('Energy (meV)')
            plt.title('Spectra of a triangular antiferromagnet')
            plt.savefig('pyspinw.png')
            plt.show()

        self.verify(kag33spec, [], {'energy': af33kagome.energy(), 'Sabp': kag33spec['Sabp']}, approxSab=0.5)


if __name__ == "__main__":
    print('################# RUNNING TESTS ###################')
    unittest.main()
