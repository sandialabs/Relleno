import os
from setuptools import setup
from setuptools.command.install import install

def readfile(filename):
    with open(filename) as f:
        return f.read()

preparse_files = ['poblano/poblano']

class SageInstall(install):
    def preparse_sage(self, prefix):
        sagename = prefix + '.sage'
        pyname = prefix + '.py'
        print('- preparsing ' + sagename)
        os.system('sage --preparse ' + sagename)
        os.system('mv ' + sagename + '.py ' + pyname)
    
    def run(self):
        print('\n\nPreparsing...')
        for prefix in preparse_files:
            self.preparse_sage(prefix)
        print('Finished preparsing!\n\n')
        install.run(self)

setup(
    name = 'poblano',
    version = readfile('VERSION'),
    description='A SageMath library for jump expansions and Jacobian matrices of conservation laws and entropy stability analysis',
    long_description = readfile('README.md'),
    license='BSD-2-Clause',
    author='Michael A. Hansen',
    author_email='mahanse@sandia.gov',
    packages = ['poblano'],
    cmdclass = {'install': SageInstall},
    url = 'https://github.com/michael-a-hansen/poblano.git'
)