#!/bin/env python
######################################################################
# Tool for querying the PhEDEx Data Service
#
# Author: Derek Feichtinger <derek.feichtinger@psi.ch>
#
# Version info: $Id: pds_simple_query.py 2256 2010-03-19 08:26:08Z dfeich $
######################################################################


import sys
import re
from optparse import OptionParser
import urllib
import xml.dom.minidom


def xmlNodeInfo(n):
    """Prints some information about the passed DOM node object"""
    print "### NODE INFO ###"
    print "# XML: %s" % (n.toprettyxml(),)
    print "# Node Type: %d" % (n.nodeType,)
    if n.nodeType == xml.dom.Node.ELEMENT_NODE:
        print "# Element Tag Name: %s" % (n.tagName,)
    print "# Node Value: %s" % (n.nodeValue,)

    if n.hasChildNodes():
        print "# Child nodes (num=%d):" % (n.childNodes.length),
        for cn in n.childNodes:
            print " %s" % (cn.nodeName),
        print

    if n.hasAttributes():
        print "# Attributes (num=%d):" % (n.attributes.length,),
        for i in range(n.attributes.length):
            print " %s," % (n.attributes.item(i-1).name),
        print "\n"


def pdsCall(call="blockreplicas",
            args=[('node',"T3_CH_PSI")],
            instance='prod',
            baseurl="http://cmsweb.cern.ch/phedex/datasvc/xml"):
    """Runs a query on the phedex data service and passes back a Document object"""
    myparam = urllib.urlencode(args)
    myurl = baseurl + '/' + instance + '/' + call + "?%s" % (myparam,)

    inf = urllib.urlopen(myurl)
    dom = xml.dom.minidom.parse(inf)
    inf.close()

    qerrors = dom.getElementsByTagName("error")

    if qerrors.length != 0:
        #xmlNodeInfo(qerrors.item(0))
        # note: not so clean. but data is in the cdata-section of the only child of the error element
        raise IOError(qerrors.item(0).firstChild.nodeValue)
        sys.exit(1)

    return dom


def recurseFetchAttr(parent,keychain,depth=0):
    # get the first key element from the keychain
    # an element looks like "name(attr1,attr2,..)"
    keychain_rgx=re.compile('^([^:]*):(.*)')
    res=keychain_rgx.match(keychain)
    if res:
        keyelem=res.group(1)
        keychain=res.group(2)
    else:
        keyelem=keychain
        keychain=''
    
    if keyelem=='':
        return
    
    # parse the key element into key and an attribute list
    keyattr_rgx=re.compile('([\w\-_]*)\(([\w\-_,]*)\)')
    res=keyattr_rgx.match(keyelem)
    if res:
        key=res.group(1)
        retattr_string=res.group(2)
    else:
        raise SyntaxError('wrong Syntax of key expression "%s"' % (keyelem))

    if retattr_string=="" or retattr_string=="ALL":
        retattr=[]
    else:
        retattr=retattr_string.split(',')

    indent=depth*'   '
    for elem in parent.getElementsByTagName(key):
        if flag_debug:
            xmlNodeInfo(elem)
        if retattr==[]:
            for i in range(elem.attributes.length):
                retattr.append(elem.attributes.item(i-1).name)
            retattr.sort()

        if len(retattr)>0:
            for attr in retattr:
                sys.stdout.write(indent)
                if flag_label:
                    sys.stdout.write("%s=" % (attr))
                sys.stdout.write("%s " % (elem.getAttribute(attr)))
            sys.stdout.write("\n")
        if keychain!='':
            recurseFetchAttr(elem,keychain,depth+1)

################################################
# MAIN


call = "blockreplicas"
qargs_string='node=T3_CH_PSI'
instance = 'prod'
retattr_string=""

usage="""usage %prog [options]
   example:
      %prog -c blockreplicas -q node=T3_CH_PSI -k 'block(name,files):replica(node)'
      %prog -c blockreplicas -q 'node=T2_CH_CSCS,block=/QCD_Pt15/Summer09-MC_31X_V3_*' -k 'block(name)'
      %prog -c filereplicas -q 'node=T2_CH_CSCS,block=/QCD_Pt15/Summer09-MC_31X_V3_7*' -k 'file(name)'
      %prog -c transferrequests -q limit=5,node=T3_CH_PSI  -k 'request(ALL):destinations(ALL):node(name):decided_by(dn)' -l
      %prog -c lfn2pfn -q node=T2_CH_CSCS,protocol=dcap,lfn=/store/mc

   Documentation: Look at http://cmsweb.cern.ch/phedex/datasvc/doc to learn about the data service calls
                  and arguments
"""
parser = OptionParser(usage=usage)
parser.add_option('-c',
                  action='store',
                  dest='call',
                  help='call/method to invoke on the data service (default: %s)' % call,
                  default=call)
parser.add_option('-d',
                  action='store_true',
                  dest='flag_debug',
                  help='debug mode (will also print full XML answer)',
                  default=False)
parser.add_option('-i',
                  action='store',
                  dest='instance',
                  help='data base instance to use (default: %s)' % instance,
                  default=instance)
parser.add_option('-k',
                  action='store',
                  dest='keychain',
                  help='elements/attributes to retrieve. Attributes to be printed are added in parentheses, e.g. -k "block(name,bytes):replica(node)". To get all attributes you can use ALL, e.g. -k "block(ALL)"')
parser.add_option('-l',
                  action='store_true',
                  dest='flag_label',
                  help='label the output fields with field names',
                  default=False)
parser.add_option('-q',
                  action='store',
                  dest='qargs_string',
                  help='colon separated list of query argument/value pairs, e.g. -q "node=ABC,block=xxx"',
                  default=qargs_string)

(options, args) = parser.parse_args()

# default options need to be set according to the specific call
default_opts={'blockreplicas' : {'keychain' : 'block(name,files,bytes)'},
              'filereplicas' : {'keychain' : 'file(name)'},
              'groups' : {'keychain' : 'group(name)'},
              'lfn2pfn' : {'keychain' : 'mapping(pfn)'},
              'links' : {'keychain' : 'link(from,to,status)'},
              'groupusage' : {'keychain' : 'node(name):group(name,node_files,node_bytes)'},
              'bounce' : {'keychain' : 'bounce(ALL)'}
              }


call=options.call


if options.keychain==None:
    if call in default_opts:
        keychain=default_opts[call]['keychain']
    else:
        keychain="undefined"
else:
    keychain=options.keychain

flag_debug=options.flag_debug
flag_label=options.flag_label
qargs_string=options.qargs_string

qargs=[]
map_rgx=re.compile('([^=]*)=([^=]*)')
for s in qargs_string.split(','):
    res=map_rgx.match(s)
    if res:
        qargs.append((res.group(1),res.group(2)))
    else:
        sys.stderr.write("Error: query argument cannot be parsed: %s\n" % s)
        sys.exit(1)

if flag_debug:
    print "call=%s keychain=%s" % (call,keychain)
    print "query arguments: %s" % (qargs,)


try:
    dom = pdsCall(call,qargs,instance)
except IOError,msg:
    sys.stderr.write("Error: %s\n" % msg)
    sys.exit(1)
except:
     sys.stderr.write("Error: Unknown error in pdsCall: %s : %s\n" % (sys.exc_info()[0],sys.exc_info()[1]))
     sys.exit(1)

if flag_debug:
    xmlNodeInfo(dom)
    #print "#################\n" + dom.toprettyxml() + "\n######################"

try:
    recurseFetchAttr(dom,keychain,0)
except SyntaxError,msg:
    sys.stderr.write("Error: %s\n" % msg)
    sys.exit(1)
except:
    sys.stderr.write("Error: Unknown error in pdsCall: %s : %s\n" % (sys.exc_info()[0],sys.exc_info()[1]))
    sys.exit(1)

dom.unlink()
sys.exit(0)
