from operator import methodcaller
import os.path

class FieldConverter:

	def __init__(self):
		self.sessionFile = ""
		self.inputFile = ""
		self.outputFile = ""

		self.processModuleList = []

		self.moduleFactory = "" # instantiate a new NekFactory for modules
		self.availiableModuleList = [] # get a list of available modules from factory as a list of names
		self.moduleList = []

		self.fieldSharedPtr = "" # FieldSharedPtr
		self.variableMap = ""

	def addSessionFile(self, sessionFile):
		assert is_file(sessionFile)
		self.sessionFile = sessionFile

	def addInputFile(self, inputFile):
		assert is_file(inputFile)
		self.inputFile = inputFile

	def addOutputFile(self, outputFile):
		self.outputFile = outputFile

	def addModule(self, moduleName):
		assert moduleName in self.availiableModuleList
		self.moduleList.append(moduleName)

	def _checkFiles(self):
		assert self.sessionFile is not None
		assert self.inputFile is not None
		assert self.outputFile is not None

	def _getFileExtenstion(self, fileName):
		assert '.' in fileName
		extenstion = fileName.split('.')[-1]
		return extenstion

	def _createModuleInstance(self, moduleType, moduleName):
		module = "" # instantiate new module
		module.SetDefaults()
		self.moduleList.append(module)

	def _createProcessModuleInstances(self):
		for processModule in self.processModuleList:
			self._createModuleInstance(eProcessModule, processModule)

	def _setModulePriorityList(self):
		self.moduleList = sorted(self.moduleList, key=methodcaller('GetModulePriority'))

	def _constructVariableMap(self):
		# recreate variable map to input into Process()

	def _runModules(self):
		for module in self.moduleList:
			module.Process(self.variableMap)

	def convert(self):
		self._checkFiles()

		self._createModuleInstance(eInputModule, _getFileExtenstion(self.sessionFile))
		self._createModuleInstance(eInputModule, _getFileExtenstion(self.inputFile))
		self._createModuleInstance(eOutputModule, _getFileExtenstion(self.outputFile))

		self._createProcessModuleInstances()

		self._setModulePriorityList()

		self._constructVariableMap()

		self._runModules()





def main():

    SessionFile = "chan3DH1D.xml"
    InputFile = "chan3DH1D.fld"
    OutputFile = "chan3DH1D.dat"

    converter = FieldConverter()

    converter.addSessionFile(SessionFile)
    converter.addInputFile(InputFile)
    converter.addOutputFile(OutputFile)

    converter.convert()
    

if __name__ == "__main__":
    main()