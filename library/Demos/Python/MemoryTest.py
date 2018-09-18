from NekPy.LibUtilities import SessionReader
from NekPy.SpatialDomains import MeshGraph
from NekPy.MultiRegions import ExpList2D
import sys
import ctypes
import gc
# import guppy

def get_refcount(coords_address):
	gc.collect()
	# return len(gc.get_referrers(ctypes.cast(coords_address, ctypes.py_object).value))
	return ctypes.c_long.from_address(coords_address).value

# def namestr(obj, namespace):
#     return [name for name in namespace if namespace[name] is obj]

def main(): 
	# hp = guppy.hpy()
	session_name = ["MemoryTest.py", "newsquare_2x2.xml"]

	session = SessionReader.CreateInstance(session_name)
	graph = MeshGraph.Read(session)
	exp = ExpList2D(session, graph)

	print("Loaded session: %s" % session.GetSessionName())
	print("Loaded MeshGraph of dimension: %d" % graph.GetMeshDimension())
	print("\n")

	coords = exp.GetCoords()
	coords1 = coords[1]
	coords = coords[0]

	# print coords

	coords_address = id(coords)
	coords_data_address = coords.ctypes.data
	coords_data_size = coords.nbytes
	print("Got coordinates...")
	print("Reference count for expansion coordinates: %d" % get_refcount(coords_address)) #1
	print "The memory address of coordinates is: {}\n".format(hex(coords_address))
	print "The memory address of data is: {}\n".format(hex(coords_data_address))
	print "The size of the data is: {}\n".format(coords_data_size)
	# print hp.iso(coords).referrers

	print("Setting PhysArray (exp.SetPhysArray())...")
	exp.SetPhysArray(coords) 
	print("exp.SetPhysArray() completed...")
	print("Reference count for expansion coordinates: %d\n" % get_refcount(coords_address)) #2
	# print hp.iso(coords).referrers
	print exp.GetPhysAddress()
	print exp.GetPhys()

	del coords
	gc.collect()
	print("del coords completed...")
	print("Reference count for expansion coordinates: %d" % get_refcount(coords_address)) #1

	print "Calculating the integral..."
	print exp.PhysIntegral()	

	print "Clearing 0x%08x size %i bytes" % (coords_data_address, coords_data_size)
	ctypes.memset(coords_data_address, 0, coords_data_size)

	print "Calculating the integral..."
	print exp.PhysIntegral()	
	

	# print "Memory has been cleared. Attempting to retrieve the data..."
	# print ctypes.cast(coords_data_address, ctypes.py_object).value

	# exp.SetPhysArray(coords1)
	# exp.DelPhys()
	# del exp
	# print("coords substituted")
	# print("Reference count for expansion coordinates: %d" % get_refcount(coords_address)) # 0


	

if __name__ == '__main__':
    main()
