# Relleno - a SageMath library for automatic jump expansion and convenient calculus                        
# Copyright 2019 National Technology & Engineering Solutions of Sandia, LLC (NTESS).                       
#                                                                                                        
# This program is free software: you can redistribute it and/or modify                                     
# it under the terms of the GNU General Public License as published by                                     
# the Free Software Foundation, either version 3 of the License, or                                        
# (at your option) any later version.                                                                      
#                                                                                                          
# This program is distributed in the hope that it will be useful,                                          
# but without any warranty; without even the implied warranty of                                           
# merchantability or fitness for a particular purpose.  See the                                            
# GNU General Public License for more details.                                                             
#                                                                                                          
# You should have received a copy of the GNU General Public License                                        
# along with this program.  If not, see <http://www.gnu.org/licenses/>.                                    
#                                                                                                          
# Questions? Contact Mike Hansen (mahanse@sandia.gov)  

import os
from setuptools import setup
from setuptools.command.install import install

def readfile(filename):
    with open(filename) as f:
        return f.read()

preparse_files = ['relleno/relleno']

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
    packages = ['relleno'],
    cmdclass = {'install': SageInstall},
    url = 'https://github.com/sandialabs/relleno.git'
)