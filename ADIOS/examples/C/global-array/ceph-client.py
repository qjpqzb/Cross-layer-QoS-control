import rados

try:
        cluster = rados.Rados(conffile='')
except TypeError as e:
        print 'Argument validation error: ', e
        raise e

print "Created cluster handle."

try:
        cluster.connect()
except Exception as e:
        print "connection error: ", e
        raise e
finally:
        print "Connected to the cluster."

print "\n\nI/O Context and Object Operations"
print "================================="

print "\nCreating a context for the 'ssd_pool' pool"
if not cluster.pool_exists('ssd_pool'):
        raise RuntimeError('No data pool exists')
#ioctx = cluster.open_ioctx('data')
ioctx = cluster.open_ioctx('ssd_pool')

#print "\nWriting object 'hw' with contents 'Hello World!' to pool 'data'."
#ioctx.write("hw", "Hello World!")
#print "Writing XATTR 'lang' with value 'en_US' to object 'hw'"
#ioctx.set_xattr("hw", "lang", "en_US")

#dou=[102.44,125.33,456.23,453.2]
#print "\nWriting object 'bm' with contents 'Bonjour tout le monde!' to pool 'data'."
#ioctx.write("bm", dou)
#print "Writing XATTR 'lang' with value 'fr_FR' to object 'bm'"
#ioctx.set_xattr("bm", "lang", "fr_FR")

print "\nContents of object 'dpot/L1'\n------------------------"
print ioctx.read("Temperature_109.913634")

#print "\n\nGetting XATTR 'lang' from object 'hw'"
#print ioctx.get_xattr("hw", "lang")

#print "\nContents of object 'bm'\n------------------------"
#print ioctx.read("bm")

#print "Getting XATTR 'lang' from object 'bm'"
#print ioctx.get_xattr("bm", "lang")


#print "\nRemoving object 'hw'"
#ioctx.remove_object("hw")

#print "Removing object 'bm'"
#ioctx.remove_object("bm")

print "\nClosing the connection."
ioctx.close()

print "Shutting down the handle."
cluster.shutdown()
