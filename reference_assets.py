"""Fetch only a public pinned reference; never send candidate sequences."""
from pathlib import Path
import hashlib,json,os,urllib.request

def ensure_reference(root):
    root=Path(root);spec=json.loads((root/'reference_manifest.json').read_text())
    cache=root/'.cache/reference/antibacterial.fasta'
    override=os.environ.get('AMP_REFERENCE_FASTA')
    paths=([Path(override)] if override else [])+[cache]+[root/p for p in spec['destinations']]
    payload=None;source=None
    for path in paths:
        if path.exists():
            content=path.read_bytes()
            if hashlib.sha256(content).hexdigest()!=spec['sha256']:
                raise ValueError('Reference hash mismatch; refusing changed cached/supplied reference')
            payload=content;source='supplied_local_file' if override and path==Path(override) else 'verified_local_cache';break
    if payload is None:
        print('Obtaining the pinned public organizer reference; no candidates are transmitted.',flush=True)
        request=urllib.request.Request(spec['url'],headers={'User-Agent':'amp-challenge-full-tier/1.1'})
        # One retry for a transient transport failure. Content/hash failure never falls back.
        for attempt in range(2):
            try:
                with urllib.request.urlopen(request,timeout=60) as response:payload=response.read(spec['bytes']+1)
                break
            except OSError:
                if attempt:raise
        source='pinned_public_https'
    if len(payload)!=spec['bytes'] or hashlib.sha256(payload).hexdigest()!=spec['sha256']:
        raise ValueError('Reference download/length/hash mismatch')
    for path in [cache]+[root/p for p in spec['destinations']]:
        path.parent.mkdir(parents=True,exist_ok=True)
        if not path.exists():path.write_bytes(payload)
        elif path.read_bytes()!=payload:raise ValueError('Conflicting source reference asset')
    return {'source':source,'sha256':spec['sha256'],'public_url':spec['url'],'candidate_data_transmitted':False}
