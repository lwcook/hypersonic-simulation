from setuptools import setup

setup(name='hypersonicvehicle',
        version='0.0.1',
        url='https://www-edc.eng.cam.ac.uk/aerotools/hypersonicvehicle',
        download_url='https://github.com/lwcook/hypersonic-vehicle/archive/0.0.1.tar.gz',
        author='Laurence W. Cook',
        author_email='lwc24@cam.ac.uk',
        install_requires=['numpy >= 1.12.1'],
        license='MIT',
        packages=['horsetailmatching'],
        test_suite='nose.collector',
        tests_require=['nose'],
        include_package_data=True,
        zip_safe=False)
