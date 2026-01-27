cmake = import('cmake')
mfem_cmake_options = cmake.subproject_options()
mfem_cmake_options.add_cmake_defines({
    'MFEM_ENABLE_EXAMPLES': 'OFF',
    'MFEM_ENABLE_TESTING': 'OFF',
    'MFEM_ENABLE_MINIAPPS': 'OFF',
    'MFEM_USE_BENCMARK': 'OFF',
    'BUILD_SHARED_LIBS': 'ON',
    'CMAKE_SKIP_INSTALL_RULES': 'ON'
})

mfem_sp = cmake.subproject(
    'mfem',
    options: mfem_cmake_options)
mfem_dep = mfem_sp.dependency('mfem')
