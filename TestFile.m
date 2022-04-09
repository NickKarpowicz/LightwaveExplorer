fid = fopen('TestFile_ExtIn.dat','rb'); 
TestFile_ExtIn = fread(fid, 1048576, 'double'); 
TestFile_ExtIn = reshape(TestFile_ExtIn,[1024 512 2]); 
fclose(fid); 
fid = fopen('TestFile_ExtOut.dat','rb'); 
TestFile_ExtOut = fread(fid, 1048576, 'double'); 
TestFile_ExtOut = reshape(TestFile_ExtOut,[1024 512 2]); 
fclose(fid); 
