message('Found CXX compiler: ' + meson.get_compiler('cpp').get_id())
message('C++ standard set to: ' + get_option('cpp_std'))

ignore_unused_args = '-Wno-unused-command-line-argument'

add_global_arguments(ignore_unused_args, language: 'cpp')

subdir('CPPC')


