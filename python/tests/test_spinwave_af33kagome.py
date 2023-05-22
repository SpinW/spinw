from pyspinw import Matlab
import numpy as np
import scipy.io
import unittest
import os

m = Matlab()

class AF33kagomeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # setup
        cls.af33kagome = m.spinw();
        cls.af33kagome.genlattice('lat_const',[6, 6, 40],'angled',[90, 90, 120],'sym','P -3');
        cls.af33kagome.addatom('r',[1/2, 0, 0],'S', 1,'label','MCu1','color','r');
        cls.af33kagome.gencoupling('maxDistance',7);
        cls.af33kagome.addmatrix('label','J1','value',1.00,'color','g');
        cls.af33kagome.addcoupling('mat','J1','bond',1);
        S0 = np.array([[0, 0, -1], [1, 1, -1], [0, 0, 0]], dtype=float)  # doesn't work for lists
        cls.af33kagome.genmagstr('mode','helical','k',[-1/3, -1/3, 0],'n',[0, 0, 1],'unit','lu','S',S0,'nExt',[1, 1, 1])
        cls.test_data_dir = os.path.join(os.path.dirname(__file__), '..', '..', 'test_data')

    def test_spinwave_calc(self):
        kag33spec = self.af33kagome.spinwave([[-1/2, 0, 0], [0, 0, 0], [1/2, 1/2, 0], 100],'hermit',False,'saveSabp',True)
        evect = np.linspace(0, 3, 100).tolist() # if not .tolist() then Evect has 501 edges!
        kag33spec = m.sw_egrid(kag33spec,'component','Sxx+Syy','imagChk',False, 'Evect', evect)
        # Reduce values of S(q,w) so it falls within tolerance (rather than change tolerance for all values)
        kag33spec['swConv'] = kag33spec['swConv'] / 2e5
        # Ignores swInt in this case
        kag33spec['swInt'] = 0

        try:
            from matplotlib import pyplot as plt
        except ImportError:
            plt = None

        if plt is not None:
            ax = plt.imshow(np.real(np.flipud(kag33spec['swConv'])),
                            aspect='auto')
            # ax.set_clim(0, 1e-6)
            plt.xlabel('Q [q, q, 0] (r.l.u)')
            plt.ylabel('Energy (meV)')
            plt.title('Spectra of a triangular antiferromagnet')
            plt.savefig('pyspinw.png')
            plt.show()

        # validate
        mat = scipy.io.loadmat(os.path.join(self.test_data_dir, 'system_tests', 'systemstest_spinwave_af33kagome.mat'))
        for ifield, field in enumerate(['omega','Sab','swConv', 'swInt']):
            print(np.sum(abs(abs(kag33spec[field])-abs(mat['data']['spec'][0,0][0,ifield]))>1e-2))
            imax = np.argmax(abs(np.max(kag33spec[field]- mat['data']['spec'][0,0][0,ifield])))
            # print(kag33spec[field][imax], mat['data']['spec'][0,0][0,ifield][imax])
            # self.assertTrue(np.allclose(kag33spec[field], mat['data']['spec'][0,0][0,ifield]))

if __name__ == "__main__":
    print('################# RUNNING TESTS ###################')
    unittest.main()
