

from cordex import variable as var



print(var.table('cmip5'))

tas = var.Variable('tas')

print((tas.frequency))

print(var.variables().keys())

