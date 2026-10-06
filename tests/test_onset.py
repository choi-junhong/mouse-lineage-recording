"""Numerical checks for the leading half-maximum onset definition."""
import unittest
import numpy as np
from fatevec.analysis import find_onset_times


class OnsetTests(unittest.TestCase):
    def setUp(self):
        self.x = np.linspace(0, 1, 1001)

    def test_quadratic_crossing_is_halfway_not_at_velocity_peak(self):
        # F=x^2 gives velocity 2x: its half-maximum crossing is exactly x=0.5.
        onset, maxima, velocity = find_onset_times({0: 0.2 + 0.8*self.x**2}, self.x)
        self.assertAlmostEqual(onset[0], 0.5, places=8)
        self.assertAlmostEqual(maxima[0], 2.0, places=8)
        self.assertGreater(self.x[np.argmax(velocity[0])], onset[0])

    def test_first_crossing_precedes_later_larger_wave(self):
        rate = (0.7*np.exp(-((self.x-0.2)/0.04)**2) +
                np.exp(-((self.x-0.75)/0.05)**2))
        fraction = np.cumsum(rate);fraction -= fraction[0];fraction /= fraction[-1]
        onset, maxima, derivative = find_onset_times({0: fraction}, self.x, window_length=31)
        self.assertTrue(0.1 < onset[0] < 0.2)
        self.assertGreater(self.x[np.argmax(derivative[0])], 0.7)
        self.assertAlmostEqual(np.interp(onset[0],self.x,derivative[0]),maxima[0]/2,places=9)
        self.assertTrue(np.all(derivative[0][self.x < onset[0]] < maxima[0]/2))

    def test_root_boundary_onset(self):
        onset, _, _ = find_onset_times({0: self.x}, self.x)
        self.assertEqual(onset[0], 0.0)

    def test_flat_or_already_pure_curves_are_undefined(self):
        onset, _, _ = find_onset_times({0: np.full(1001,0.3),1: np.ones(1001)},self.x)
        self.assertTrue(np.isnan(onset[0]))
        self.assertTrue(np.isnan(onset[1]))

    def test_root_baseline_does_not_change_normalized_velocity(self):
        onset, _, velocity = find_onset_times({0:self.x**2,1:0.8+0.2*self.x**2},self.x)
        self.assertAlmostEqual(onset[0],onset[1],places=9)
        np.testing.assert_allclose(velocity[0],velocity[1],atol=1e-9)

    def test_grid_units_scale_onset_and_velocity_consistently(self):
        a,m,_=find_onset_times({0:self.x**2},self.x)
        b,n,_=find_onset_times({0:self.x**2},7*self.x)
        self.assertAlmostEqual(b[0],7*a[0],places=9)
        self.assertAlmostEqual(n[0],m[0]/7,places=9)

    def test_invalid_grids_and_windows_raise(self):
        with self.assertRaises(ValueError):
            find_onset_times({0:self.x},self.x,window_length=300)
        with self.assertRaises(ValueError):
            find_onset_times({0:self.x},self.x**2)
        with self.assertRaises(ValueError):
            find_onset_times({0:self.x},self.x,window_length=1003)


if __name__ == '__main__':
    unittest.main()
