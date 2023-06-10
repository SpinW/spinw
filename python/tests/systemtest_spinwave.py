from pyspinw import Matlab
from libpymcr.MatlabProxyObject import wrap
import numpy as np
import scipy.io
import unittest
import copy
import os

m = Matlab()

class SystemTest_Spinwave(unittest.TestCase):
    """
    Port to Python of Matlab `sw_tests.system_tests.systemtest_spinwave` classes
    without code to generate reference files (will use Matlab generated files)
    """

    @classmethod
    def setUpClass(cls):
        try:
            import matplotlib.pyplot
            cls.plt = matplotlib.pyplot
        except ImportError:
            cls.plt = None
        curdir = os.path.dirname(__file__)
        cls.ref_data_dir = os.path.abspath(os.path.join(curdir, '..', '..', 'test_data', 'system_tests'))
        cls.relToll = 0.01
        cls.absToll = 1.0e-6
        cls.bigVal = 1.0e8
        cls.tolSab = 0.05
        m.swpref().fid = 0   # Suppress printing to make test output less verbose


    @classmethod
    def get_hash(cls, obj):
        # Uses Matlab's (undocumented) hash to be consistent with generated reference data
        engine = wrap(m._interface.call('java.security.MessageDigest.getInstance', ['MD5']), m._interface)
        engine.update(m.getByteStreamFromArray(obj))
        return m.typecast(engine.digest(), 'uint8')


    @classmethod
    def harmonize(cls, inp):
        # Harmonizes saved mat file and new output data
        # Removes singleton dimensions and nested single structures, converts some types to match
        if isinstance(inp, np.ndarray):
            if len(inp.dtype) > 1 and all([inp.dtype[ii]=='O' for ii in range(len(inp.dtype))]):
                # Convert from dtype array to dict if array is all objects
                return {ky:cls.harmonize(inp[ky]) for ky in inp.dtype.names}
            if len(inp.shape) > 0 and inp.size == 1:
                return cls.harmonize(inp[0])
            if len(inp.shape) > 1 and any([s==1 for s in inp.shape]):
                return cls.harmonize(np.squeeze(inp))
            if inp.dtype == 'O' and inp.size > 1:
                # Convert object array to list
                return [inp[ii] for ii in range(inp.size)]
        elif isinstance(inp, dict):
            return {ky:cls.harmonize(vl) for ky, vl in inp.items()}
        elif isinstance(inp, np.str_):
            return str(inp)
        elif isinstance(inp, np.uint8):
            return np.float64(inp)
        elif isinstance(inp, np.int16):
            return np.float64(inp)
        return inp


    def load_ref_dat(self, filename):
        self.reference_data = self.harmonize(scipy.io.loadmat(os.path.join(self.ref_data_dir, filename)))


    def get_fieldname(self, pars):
        if not pars:
            pars = 'data'
        if isinstance(pars, str):
            return pars
        return f'd{m.dec2hex(self.get_hash(pars))}'


    def sanitize(self, indat):
         out = copy.deepcopy(indat)
         out[np.where(np.abs(out) > self.bigVal)] = 0.0
         return out


    def verifyIsEqual(self, test_data, ref_data, key=''):
        # Equivalent to Matlab recursive `verifyThat(a, IsEqual(b))` with absolute and relative bounds
        self.assertIs(type(test_data), type(ref_data), msg=f'Item {key} type mismatch')
        if isinstance(test_data, dict):
            self.assertEqual(test_data.keys(), ref_data.keys())
            for ky in test_data.keys():
                self.verifyIsEqual(test_data[ky], ref_data[ky], ky if not key else f'{key}.{ky}')
        elif isinstance(test_data, list):
            for ii in range(len(test_data)):
                self.verifyIsEqual(test_data[ii], ref_data[ii], f'[{ii}]' if not key else f'{key}[{ii}]')
        elif isinstance(test_data, np.ndarray):
            np.testing.assert_allclose(self.sanitize(test_data), self.sanitize(ref_data),
                                       rtol=self.relToll, atol=self.absToll, equal_nan=True,
                                       err_msg=f'Item {key} are not close', verbose=True)
        else:
            self.assertEqual(test_data, ref_data)


    def approxMatrix(self, actual, expected, frac_not_match):
        # Checks if two arrays are approximately the same with most entries equal but a fraction not
        if isinstance(actual, list):
            return [self.approxMatrix(xx, yy, frac_not_match) for xx, yy in zip(actual, expected)]
        diff = np.abs(actual - expected)
        rel_diff = np.divide(diff, expected, out=np.zeros_like(expected), where=expected!=0)
        # Possible alternative:
        #return np.where((diff > self.absToll) & (rel_diff > self.relToll), expected, actual)
        frac = np.where((diff > self.absToll) & (rel_diff > self.relToll))[0].size / diff.size
        return expected if frac < frac_not_match else actual


    def verify_eigval_sort(self, actual, expected, nested=0):
        # Checks if eigenvalues match, if not try different sort permutations
        if isinstance(actual, list):
            vv = (self.verify_eigval_sort(xx, yy, nested) for xx, yy in zip(actual, expected))
            return (sum(zz, start=zz[0]) for zz in zip(*vv))
        if not np.allclose(actual, expected, rtol=self.relToll, atol=self.absToll, equal_nan=True):
            sort_ax = np.where(np.array(actual.shape) > 1)[0][0]
            if nested > 1:
                actual = np.sort(np.abs(actual), axis=sort_ax)
                expected = np.sort(np.abs(expected), axis=sort_ax)
            else:
                actual = np.sort(np.real(actual), axis=sort_ax)
                expected = np.sort(np.real(expected), axis=sort_ax)
            if nested < 2:
                actual, expected = self.verify_eigval_sort(actual, expected, nested + 1)
        return actual, expected


    def verify(self, spec, pars, extrafields=None, approxSab=False):
        # This is the Matlab `generate_or_verify` method without the reference data generation code
        ref_data = self.reference_data[self.get_fieldname(pars)]
        # There are some type-mismatch (strings/np.str_ and ints/floats) in the input data, ignore for now
        ref_data.pop('input')
        #test_data = {'input': m.struct(self.swobj)}
        test_data = {}
        omega, ref_data['spec'][0] = self.verify_eigval_sort(spec['omega'], ref_data['spec'][0])
        test_data['spec'] = [omega, spec['Sab']]
        if 'swConv' in spec:
            test_data['spec'].append(spec['swConv'])
            if 'swInt' in spec and spec['swInt']:
            # Matlab code is not explicit that if 'swInt' in included, 'swConv' is also, but it is implicit
                test_data['spec'].append(spec['swInt'])
        else:
            assert 'swInt' not in spec
        if extrafields is not None:
            test_data.update(extrafields)
        if approxSab:
            tolSab = approxSab if isinstance(approxSab, float) else self.tolSab
            test_data['spec'][1] = self.approxMatrix(spec['Sab'], ref_data['spec'][1], tolSab)
            if len(test_data) == 4:
                test_data['spec'][3] = self.approxMatrix(spec['swInt'], ref_data['spec'][3], tolSab)
            if 'Sabp' in test_data:
                test_data['Sabp'] = self.approxMatrix(test_data['Sabp'], ref_data['Sabp'], tolSab)
            if 'V' in test_data:
                test_data['V'] = self.approxMatrix(test_data['V'], ref_data['V'], tolSab)
        self.verifyIsEqual(self.harmonize(test_data), ref_data)


    def verify_generic(self, data, fieldname):
        return self.verifyIsEqual(data, self.reference_data[fieldname])
