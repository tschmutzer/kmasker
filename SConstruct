#!python
env = Environment()
opts = Variables()
opts.Add(PathVariable('PREFIX', 'Directory to install under', '/usr', PathVariable.PathIsDir))
opts.Update(env)
Help(opts.GenerateHelpText(env))

env.Default('install')

# Here are our installation paths:
idir_prefix = '$PREFIX'
#idir_lib    = '$PREFIX/lib'
#idir_bin    = '$PREFIX/bin'
#idir_etc    = '$PREFIX/etc'
#idir_data   = '$PREFIX/share'
Export('env idir_prefix')

SConscript('src/Sconscript', exports=['env', 'opts'])
env.Install(idir_prefix, 'bin')
env.Install(idir_prefix, 'lib')
env.Install(idir_prefix, 'etc')
env.Install(idir_prefix, 'share')

env.Alias('install', ['build', idir_prefix])