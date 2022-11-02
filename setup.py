import setuptools

def readme():
    try:
        with open('readme.md') as f:
            return f.read()
    except:
        pass

setuptools.setup(
    name='MultiGapFilling',
    version='1.0.3',
	license='GPLv3+',
    author="Lucas-Moffat Antje",
    author_email='Antje.Lucas-Moffat@dwd.de',
    packages=setuptools.find_packages(),
    include_package_data=True,    
    description='Multiple gap-filling tool for or eddy covariance datasets',
    long_description=readme(),
    long_description_content_type='text/markdown',
    url='https://www.icos-cp.eu/',
    project_urls={
            'Source':'https://github.com/claudiodonofrio/mgf',
			'Documentation':'https://doi.org/10.1016/j.agrformet.2022.109114'},
    install_requires=['matplotlib','scipy'],
    classifiers=[
        'Programming Language :: Python :: 3',
		'Development Status :: 5 - Production/Stable',		
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Intended Audience :: Science/Research',
        'Operating System :: OS Independent',
        'Natural Language :: English',
        'Topic :: Scientific/Engineering :: Atmospheric Science',
		'Topic :: Scientific/Engineering :: Artificial Intelligence',
		'Topic :: Scientific/Engineering :: Information Analysis',
        'Topic :: Utilities',
    ],
)






