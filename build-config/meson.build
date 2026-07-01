subdir('mfem')
subdir('libconfig')
subdir('CLI11')
subdir('magic_enum')

if get_option('build_python')
    subdir('python')
    subdir('pybind')
endif

if get_option('build_python')
    stroid_pkg_dir    = py_installation.get_install_dir() / 'stroid'
    stroid_includedir = stroid_pkg_dir / 'include'
    stroid_libdir     = stroid_pkg_dir / 'lib'
    stroid_pcdir      = stroid_libdir / 'pkgconfig'
else
    stroid_includedir = get_option('includedir')
    stroid_libdir     = get_option('libdir')
    stroid_pcdir      = get_option('libdir') / 'pkgconfig'
endif