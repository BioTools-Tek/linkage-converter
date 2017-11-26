#!/usr/bin/env python

from distutils.core import setup

setup(name = 'Linkage2Allegro',
      version = '2017.3',
      description = 'Converts between linkage output formats {Genehunter, Merlin, Simwalk, Swiftlink} --> Allegro ',
      url = 'https://github.com/BioTools-Tek/linkage-converter',
      author = 'Mehmet Tekman',
      author_email = 'mtekman89@gmail.com',
      license = 'GPL-3',
      keywords = 'linkage haplotype genehunter allegro simwalk swiftlink merlin conversion',
      package_dir = {'linkage2allegro' : 'src'},
      packages = ['linkage2allegro'],
      scripts = ['bin/linkage2allegro']
)
