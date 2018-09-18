import subprocess
import os.path

class Converter:

	def __init__(self):
		# to be obtained from cmake configuration
		self.FieldConvertPath = "/home/emi/FYP/nektar++/build/utilities/FieldConvert/FieldConvert"
		self.ModuleList = []
		self.SessionFile = ""
		self.InputFile = ""
		self.OutputFile = ""
		self.OutputFileOverride = False

	def addSessionFile(self, SessionFile):
		self.SessionFile = SessionFile

	def addInputFile(self, InputFile):
		self.InputFile = InputFile

	def addOutputFile(self, OutputFile, override=False):
		self.OutputFileOverride = override
		if os.path.isfile(OutputFile) and self.OutputFileOverride is False:
			print("The output file already exists")
		else:
			self.OutputFile = OutputFile

	def convert(self):
		p = subprocess.Popen([self.FieldConvertPath, self.SessionFile, self.InputFile, self.OutputFile], stdout=subprocess.PIPE)
		print p.communicate()


def main():

    SessionFile = "chan3DH1D.xml"
    InputFile = "chan3DH1D.fld"
    OutputFile = "chan3DH1D.dat"

    converter = Converter()

    converter.addSessionFile(SessionFile)
    converter.addInputFile(InputFile)
    converter.addOutputFile(OutputFile)

    converter.convert()
    

if __name__ == "__main__":
    main()