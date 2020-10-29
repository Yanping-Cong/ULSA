from setuptools import setup, find_packages
import os
REQUIRES = ['numpy', 'scipy', 'matplotlib', 'h5py', 'caput', 'healpy', 'astropy']
setup(name='ULSA',
      version='0.1',
      description='The Ultral-Long wavelength Sky model with Absorption',
      #long_description=readme(),
      install_requires=REQUIRES,
      classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 2.7',
        'Topic :: Text Processing :: Linguistic',
      ],
      keywords='The ULSA',
      url='https://github.com/Yanping-Cong/ULSA',
      author='Yanping Cong',
      author_email='ypcong@bao.ac.cn',
      license='MIT',
      #packages=['funniest/funniest'],
      packages = find_packages(),
      #install_requires=[
      #    'markdown',
      #],
      dependency_links = [
          'https://github.com/zuoshifan/caput/tree/zuo/develop'
          ],
      include_package_data=True,
      zip_safe=False)
