#!/usr/bin/env python

import os

try: from setuptools import setup
except: from distutils.core import setup

setup(name = 'actin',
      version = '1.2.3',
      description = 'Activity Indices Calculator',
      url = 'http://github.com/gomesdasilva/actin',
      download_url = 'https://github.com/gomesdasilva/ACTIN/archive/v1.0.tar.gz',
      author = 'Joao Gomes da Silva',
      author_email = 'Joao.Silva@astro.up.pt',
      license = 'MIT',
      keywords = ['astronomy', 'activity', 'fits', 'harps', 'harps-n'],
      packages = ['actin'],
      entry_points = {
        "console_scripts": ['actin = actin.actin:main']
        },
      include_package_data = True,
      install_requires = ['appdirs']
      )


# This runs ACTIN and gives location of config file
os.system('actin -cfg True')

#import appdirs
#cfg_dir = appdirs.user_config_dir('ACTIN')
#os.chmod(cfg_dir, 0o777) # make directory read/write by everyone
