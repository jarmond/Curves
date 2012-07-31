from distutils.core import setup

setup(name = 'Curves',
      version = '0.4.0',
      description = 'Curves Python Module',
      author = 'Jonathan Armond',
      author_email = 'j.w.armond@warwick.ac.uk',
      packages = ['curves'],
      scripts = ['scripts/calcdsens', 'scripts/cconvert', 'scripts/cprocess',
                 'scripts/peakfind', 'scripts/stepfind', 'scripts/align',
                 'scripts/wavelet', 'scripts/peak_analyse']
     )

