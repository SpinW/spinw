__author__ = "github.com/wardsimon"
__version__ = "0.0.1"

import sys
import requests

org = "spinw"
repo = "spinw"
if len(sys.argv) > 1:
    repo = sys.argv[1]

releases = requests.get(f"https://api.github.com/repos/{org}/{repo}/releases").json()

header = f"<!DOCTYPE html>\n<html>\n<head>\n<title>Links for {repo} (alpha)</title>\n</head>\n<body>\n<h1>Links for {repo}</h1>"
body = ""
for release in releases:
    asset_url = release["assets_url"]
    assets = requests.get(asset_url).json()
    for asset in assets:
        if asset["name"].endswith(".whl"):
            name = asset["name"][:-4]
            url = asset["browser_download_url"]
            body += f'<a href="{url}" data-requires-python="&gt=3.7,&lt4.0">{name}</a><br/>\n'
footer = "</body>\n</html>"

content = "\n".join([header, body, footer])
with open("index.html", "w") as fid:
    fid.write(content)
