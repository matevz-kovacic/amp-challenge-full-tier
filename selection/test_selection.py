"""Regression tests for the packaged rank helper and frozen assembly algorithms."""
import unittest
import numpy as np
from replay_selection import rankdata,borda,mmr,apportion,RFLev

class SelectionTests(unittest.TestCase):
    def test_average_ties_and_direction(self):
        np.testing.assert_array_equal(rankdata([9,3,3,6]),[4,1.5,1.5,3])
        specs=[dict(id='risk',family='risk',direction='minimize'),dict(id='preference',family='other',direction='maximize')]
        _,_,_,q,_,_=borda([[0,9],[1,3],[1,3]],specs)
        np.testing.assert_array_equal(q,[1,.25,.25])
    def test_missing_values_fail(self):
        with self.assertRaises(ValueError):borda([[np.nan]], [dict(id='a',family='a',direction='minimize')])
    def test_constants_and_copies(self):
        specs=[dict(id='a',family='a',direction='maximize'),dict(id='b',family='a',direction='maximize',implementation='a'),dict(id='c',family='c',direction='maximize')]
        _,_,_,q,active,omitted=borda([[1,1,4],[3,3,4]],specs)
        self.assertEqual(active,[0]);self.assertEqual(len(omitted),2);np.testing.assert_array_equal(q,[0,1])
    def test_singleton(self):
        _,_,_,q,_,_=borda([[8]],[dict(id='a',family='a',direction='maximize')])
        self.assertEqual(q.tolist(),[1.])
    def test_direction_reverse_and_palindrome(self):
        s='ACDEFGHI';pal='ACDEEDCA'
        self.assertEqual(len(set([s,s,s[::-1],pal,pal[::-1]])),3)
        self.assertLess(RFLev.normalized_similarity(s,s[::-1]),1)
        chosen,_=mmr([0,1],['ACDEFGHI','ACDEFGHK'],np.array([.9,.9]),2,.9)
        self.assertEqual(chosen,[0,1])
    def test_apportion_exact_capacity(self):
        np.testing.assert_array_equal(apportion([2,2,2],2),[1,1,0])
        np.testing.assert_array_equal(apportion([0,1,9],30),[0,1,9])

if __name__=='__main__':unittest.main()
