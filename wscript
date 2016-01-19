VERSION = '0.0.1'
APPNAME = 'betweenness'

srcdir = '.'
blddir = 'build'

def set_options(ctx):
  ctx.tool_options('compiler_cxx')

def configure(ctx):
  ctx.check_tool('compiler_cxx')
  ctx.check_cxx(lib = ['glog', 'gflags'], uselib_store = 'common')
  ctx.env.CXXFLAGS += ['-O2', '-Wall', '-g', '-I/usr/local/include/', '-std=c++0x']
  ctx.env.append_value('LINKFLAGS', ['-L/usr/local/lib/'])


def build(bld):
  bld(features     = 'cxx cprogram',
      source       = 'test.cc',
      target       = 'test',
      includes     = ['.'],
      lib          = ['gflags', 'glog'])


