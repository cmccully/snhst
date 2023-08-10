from setuptools import setup

setup(name='snhst',
      author=['Curtis McCully'],
      author_email=['cmccully@lco.global'],
      version='0.10.0',
      packages=['snhst'],
      setup_requires=[],
      install_requires=['numpy<=1.22', 'astropy', 'matplotlib', 'drizzlepac', 'astroscrappy', 'reproject',
                        'scipy', 'crds', 'sep'],
      tests_require=[],
      scripts=[],
      entry_points={'console_scripts': ['reduce_hst_data=snhst.reduce_hst_data:run']}
      )
