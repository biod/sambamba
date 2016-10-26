template = """package:
  name: sambamba
  version: '{version}'

source:
  fn: sambamba_v{version}.tar.bz2
  url: {linux_url} # [linux]
  md5: {linux_md5} # [linux]
  url: {osx_url} # [osx]
  md5: {osx_md5} # [osx]

build:
  number: 0

requirements:
  build:
  run:
    - samtools # required for mpileup
    - bcftools # required for mpileup

test:
  commands:
    - sambamba view

about:
  home: https://github.com/lomereiter/sambamba
  license: GPLv2
  summary: Tools for working with SAM/BAM data"""

import json
from urllib2 import urlopen
from hashlib import md5

latest_release = json.loads(urlopen("https://api.github.com/repos/lomereiter/sambamba/releases").read())[0]
sambamba_version = latest_release['tag_name'][1:]

downloads = {}
for asset in latest_release['assets']:
    url = asset['browser_download_url']
    platform = asset['name'].split(sambamba_version)[1].split(".")[0][1:]
    downloads[platform] = url

def md5sum(download):
    h = md5()
    h.update(urlopen(download).read())
    return h.hexdigest()

linux_md5 = md5sum(downloads['linux'])
osx_md5 = md5sum(downloads['osx'])

print template.format(version=sambamba_version, 
                linux_url=downloads['linux'], 
                osx_url=downloads['osx'],
                linux_md5=linux_md5, osx_md5=osx_md5)
