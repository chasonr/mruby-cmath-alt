# This `mruby-cmath-alt` gem uses C99 _Complex features and GCC extensions,
# but not complex.h.
# You need a version of GCC that supports C99+.
MRuby::Gem::Specification.new('mruby-cmath-alt') do |spec|
  spec.license = 'MIT'
  spec.author  = 'mruby developers'
  spec.summary = 'standard Math module with complex'
  spec.add_dependency 'mruby-complex', :core => 'mruby-complex'
end
