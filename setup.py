from setuptools import setup, find_packages

REQUIRES = [
          'numpy == 1.15.0',
          'scikit-image',
          'scikit-learn',
          'nested_lookup',
          'read-roi',
          'scipy',
          'matplotlib',
          'Pillow',
          'matplotlib-venn',
          'pandas',
          'seaborn',
]

setup(
      name='rnaloc',
      version='0.1.3',
      description='Code to analyze RNA localization in c. elegans.',
      url='https://github.com/muellerflorian/parker-rna-loc-elegans',
      author='Florian MUELLER',
      author_email='muellerf.research@gmail.com',
      license='MIT',
      packages=find_packages(),
      include_package_data=True,
      install_requires=REQUIRES,
      zip_safe=False
      )
