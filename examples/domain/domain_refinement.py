from cordex import domain as dm



eur44 = dm.domain('EUR-44')
eur22 = dm.domain('EUR-22')
eur11 = dm.domain('EUR-11')

# use of refinement funtion
eur44_ref = eur44.refine(2.0)

# compare eur44_ref with EUR-22 from the table
print(eur22 == eur44_ref)

# demonstrate simple domain math.
print(eur22 == 0.5 * eur44 )
print(eur11 == 0.25 * eur44 )
print(eur11 == 0.5 * eur22 )
print(eur44 == 4 * eur11 )

