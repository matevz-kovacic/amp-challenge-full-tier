"""Public-release regression checks: exact novelty boundary and pinned reference."""
from pathlib import Path
import hashlib,importlib.util,json,os,sys,tempfile,unittest
from unittest.mock import patch
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from reference_assets import ensure_reference
from rapidfuzz.distance import Indel,Levenshtein as UnitEdit
spec=importlib.util.spec_from_file_location('portable_checks',ROOT/'scripts/portable_sequence_checks.py')
checks=importlib.util.module_from_spec(spec);spec.loader.exec_module(checks)

class PublicReleaseTests(unittest.TestCase):
    def test_indel_boundary_is_distinct_from_unit_edit(self):
        s='ACDEFGHIKL';exact='ACDEFGHIAA';over='ACDEFGHIKA'
        self.assertEqual(Indel.normalized_similarity(s,exact),.8)
        self.assertGreater(Indel.normalized_similarity(s,over),.8)
        checks._veritfy_max_simularity({s},{exact})
        with self.assertRaises(ValueError):checks._veritfy_max_simularity({s},{over})
        # An insertion demonstrates the two normalizations are not interchangeable.
        self.assertNotEqual(Indel.normalized_similarity(s,s+'A'),UnitEdit.normalized_similarity(s,s+'A'))
    def test_reversed_peptides_remain_distinct(self):
        s='ACDEFGHI';reverse=s[::-1]
        self.assertEqual(len({s,reverse}),2)
        checks._verify_no_overlap({s},{reverse})
        checks._veritfy_max_simularity({s},{reverse})
    def test_reference_offline_copy_cache_and_hash_failure(self):
        with tempfile.TemporaryDirectory(prefix='amp-public-reference-test-') as folder:
            root=Path(folder);payload=b'>public-reference\nACDEFGHI\n';source=root/'provided.fasta';source.write_bytes(payload)
            cfg={'url':'https://invalid.example.test/no-network','bytes':len(payload),'sha256':hashlib.sha256(payload).hexdigest(),'destinations':['source/antibacterial.fasta']}
            (root/'reference_manifest.json').write_text(json.dumps(cfg))
            with patch.dict(os.environ,{'AMP_REFERENCE_FASTA':str(source)}),patch('urllib.request.urlopen') as network:
                result=ensure_reference(root);self.assertEqual(result['source'],'supplied_local_file');network.assert_not_called()
            self.assertEqual((root/'source/antibacterial.fasta').read_bytes(),payload)
            with patch.dict(os.environ,{},clear=True),patch('urllib.request.urlopen') as network:
                self.assertEqual(ensure_reference(root)['source'],'verified_local_cache');network.assert_not_called()
                (root/'.cache/reference/antibacterial.fasta').write_bytes(b'changed')
                with self.assertRaisesRegex(ValueError,'Reference hash mismatch'):ensure_reference(root)
                network.assert_not_called()

if __name__=='__main__':unittest.main()
