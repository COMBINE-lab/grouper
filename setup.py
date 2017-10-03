from setuptools import setup

setup(name='biogrouper',
      version='0.1.2',
      scripts=['bin/Grouper'],
      description='Graph-based clustering and annotation for improved de novo transcriptome analysis',
      url='https://github.com/COMBINE-lab/Grouper',
      author='Laraib Malik, Fatemeh Almodaresi, Rob Patro',
      author_email='rob.patro@cs.stonybrook.edu',
      license='BSD with attribution',
      packages=['grouper'],
      install_requires=[
          'PyYAML',
          'coloredlogs',
          'click',
          'networkx',
          'numpy',
          'pandas',
          'tqdm'
      ],
      zip_safe=False)
