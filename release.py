import argparse
import json
import os
import re
import sys
import subprocess

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--github', action='store_true', help='Release on Github')
    parser.add_argument('--create_tag', action='store_true', help='Create git tag if needed')
    parser.add_argument('--pypi', action='store_true', help='Release on PyPI')
    parser.add_argument('--notest', action='store_true', help='Actually send/upload')
    parser.add_argument('--token', action='store', help='Github token to access repo')
    parser.add_argument('--version_check', action='store_true', help='Check version strings')
    args = parser.parse_args()

    if args.version_check:
        file_ver, _ = _version_check()
        print(f'Version string "{file_ver}" in files match')

    token = args.token
    if token is None and 'GITHUB_TOKEN' in os.environ:
        token = os.environ['GITHUB_TOKEN']

    test = not args.notest
    if args.github:
        release_github(test, args.create_tag, token)

    if args.pypi:
        release_pypi(test, token)


def release_github(test=True, create_tag=False, token=None):
    rv = subprocess.run(['git', 'describe', '--tags', '--dirty', '--always', '--long'],
            capture_output=True)
    if rv.returncode != 0:
        raise Exception(f'During git describe, got this error: {rv.stderr}')
    git_ver = rv.stdout.decode().strip()
    file_ver, changelog = _version_check()
    if 'g' in git_ver and create_tag:
        # Not in a release, create a new tag
        rv = subprocess.run(['git', 'tag', file_ver], capture_output=True)
        if rv.returncode != 0:
            raise Exception(f'During tag, git returned this error: {rv.stderr}')
        git_ver = file_ver
    elif git_ver != file_ver:
        raise Exception(f'version mismatch! __version__: {git_ver}; files: {file_ver}')

    desc = re.search('# \[v[0-9\.]*\]\(http.*?\)\n(.*?)# \[v[0-9\.]*\]', changelog,
                     re.DOTALL | re.MULTILINE).groups()[0].strip()
    payload = {
        "tag_name": git_ver,
        "target_commitish": "main",
        "name": git_ver,
        "body": desc,
        "draft": False,
        "prerelease": True
    }
    if test:
        print(payload)
    else:
        upload_url = release_exists(git_ver, retval='upload_url', token=token)
        if not upload_url:
            upload_url = _create_gh_release(payload, token)
        else:
            upload_url = re.search('^(.*)\{\?', upload_url).groups()[0]
        _upload_assets(upload_url, token)


def release_pypi(test=True, token=None):
    # Downloads wheels from github and upload to PyPI
    import requests
    response = requests.get(
        'https://api.github.com/repos/spinw/spinw/releases')
    # Get the latest release
    releases = response.json()
    ids = [r['id'] for r in releases]
    latest = [r for r in releases if r['id'] == max(ids)][0]
    # Creates a custom wheelhouse folder
    try:
        os.mkdir('twine_wheelhouse')
    except FileExistsError:
        pass
    # Loops through assets and downloads all the wheels
    headers = {"Accept":"application/octet-stream"}
    for asset in latest['assets']:
        if asset['name'].endswith('whl'):
            print('Downloading %s' % (asset['name']))
            localfile = os.path.join('twine_wheelhouse', asset['name'])
            download_github(asset['url'], localfile, token)
    if not test:
        rv = subprocess.run(['twine', 'upload', 'twine_wheelhouse/*'], capture_output=True)
        if rv.returncode != 0:
            raise Exception(f'During upload, twine returned this error: {rv.stderr}')


def release_exists(tag_name, retval='upload_url', token=None):
    import requests
    headers = {}
    if token is not None:
        headers = {"Authorization": "token " + token}
    response = requests.get(
        'https://api.github.com/repos/spinw/spinw/releases',
        headers=headers)
    if response.status_code != 200:
        raise RuntimeError('Could not query Github if release exists')
    response = json.loads(response.text)
    desired_release = [v for v in response if v['tag_name'] == tag_name]
    if desired_release:
        return desired_release[0][retval]
    else:
        return False


def download_github(url, local_filename=None, token=None):
    import requests
    headers = {"Accept":"application/octet-stream"}
    if token is not None:
        headers["Authorization"] = "token " + token
    if not local_filename:
        local_filename = url.split('/')[-1]
    with requests.get(url, stream=True, headers=headers) as r:
        with open(local_filename, 'wb') as f:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)
    return local_filename


def _version_check():
    with open('CHANGELOG.md') as f:
        changelog = f.read()
    with open('CITATION.cff') as f:
        citation = f.read()
    cl_ver = re.findall('# \[(.*)\]\(http', changelog)[0]
    cit_ver = 'v' + re.findall('\nversion: "(.*)"', citation)[0]
    if cl_ver != cit_ver:
        raise Exception(f'version mismatch! CHANGELOG.md: {cl_ver}; CITATION.cff: {cit_ver}')
    return cl_ver, changelog


def _create_gh_release(payload, token):
    assert token is not None, 'Need token for this action'
    import requests
    response = requests.post(
        'https://api.github.com/repos/spinw/spinw/releases',
        data=json.dumps(payload),
        headers={"Authorization": "token " + token})
    print(response.text)
    if response.status_code != 201:
        raise RuntimeError('Could not create release')
    upload_url = re.search('^(.*)\{\?', json.loads(response.text)['upload_url']).groups()[0]
    return upload_url


def _upload_assets(upload_url, token):
    assert token is not None, 'Need token for this action'
    import requests
    for wheelpath in ['build', 'dist', 'wheelhouse']:
        wheelfile = os.path.basename(wheelpath)
        print(f'Uploading wheel {wheelpath}')
        with open(wheelpath, 'rb') as f:
            upload_response = requests.post(
                f"{upload_url}?name={wheelfile}",
                headers={"Authorization": "token " + token,
                         "Content-type": "application/octet-stream"},
                data=f.read())
            print(upload_response.text)
    return None


if __name__ == '__main__':
    main()
