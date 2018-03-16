from setuptools import setup

setup(name='biogrouper',
      python_requires="<3.0",
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
          'networkx==1.11',
          'numpy',
          'pandas',
          'tqdm'
          'statistics'
      ],
      zip_safe=False)
