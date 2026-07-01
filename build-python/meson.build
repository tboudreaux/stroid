if get_option('build_python')
    message('Building Python bindings...')

    stroid_py_deps = [
        py_dep,
        pybind11_dep,
        stroid_dep
    ]

    if host_machine.system() == 'darwin'
        stroid_ext_rpath = '@loader_path/lib'
    else
        stroid_ext_rpath = '$ORIGIN/lib'
    endif

    py_sources =  [
        meson.project_source_root() + '/src/python/bindings.cpp',
        meson.project_source_root() + '/src/python/config/bindings.cpp',
        meson.project_source_root() + '/src/python/exceptions/bindings.cpp',
        meson.project_source_root() + '/src/python/IO/bindings.cpp',
        meson.project_source_root() + '/src/python/refinement/bindings.cpp',
        meson.project_source_root() + '/src/python/utils/bindings.cpp',
    ]

    py_mod = py_installation.extension_module(
        '_stroid',
        sources: py_sources,
        dependencies: stroid_py_deps,
        install: true,
        link_args: stroid_ext_rpath_args,
        build_rpath: stroid_ext_rpath,
        install_rpath: stroid_ext_rpath,
        subdir: 'stroid',
    )

    py_installation.install_sources(
        meson.project_source_root() + '/src/python/stroid/__init__.py',
        subdir: 'stroid',
    )

else
    message('Python bindings disabled')
endif