oo = open('a', 'w')

varname = 'fTtauID%(id)s'

print >> oo, "//defining the tauid variables, by Hamed"
for idname in process.patTaus.tauIDSources.parameterNames_():
    print >> oo, ("std::auto_ptr<std::vector<float> >  " + varname + ";") % {"id":idname}
print >> oo, "//END"


print >> oo, "//tauid addProducts, by Hamed"
for idname in process.patTaus.tauIDSources.parameterNames_():
    print >> oo, ("addProduct(\"%(id)s\", typeid(*"+ varname  +"));") % {"id":idname}
print >> oo, "//END"

print >> oo, "//tauid puts, by Hamed"
for idname in process.patTaus.tauIDSources.parameterNames_():
    print >> oo, ("e.put("+varname+",fullName(\"%(id)s\"));") % {"id":idname}
print >> oo, "//END"

print >> oo, "//tauid resets, by Hamed"
for idname in process.patTaus.tauIDSources.parameterNames_():
    print >> oo, (varname+".reset(new std::vector<float>);") % {"id":idname}
print >> oo, "//END"    
  
print >> oo , "//tauid push_backs, by Hamed"
for a in process.patTaus.tauIDSources.parameterNames_():
    print >> oo, ("fTtauID%(id)s->push_back( lepton.tauID(\"%(id)s\") );") % {"id":a}
print >> oo, "//END"

oo.close()
