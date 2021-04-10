% This make.m is for MATLAB and OCTAVE under Windows, Mac, and Unix

try
	Type = ver;
	% This part is for OCTAVE
	if(strcmp(Type(1).Name, 'Octave') == 1)
		mex libsvmread.cpp
		mex libsvmwrite.cpp
		mex svmtrain.cpp ../svm.cpp svm_model_matlab.cpp
		mex svmpredict.cpp ../svm.cpp svm_model_matlab.cpp
	% This part is for MATLAB
	% Add -largeArrayDims on 64-bit machines of MATLAB
	else
		mex CFLAGS="\$CFLAGS -std=c99" -largeArrayDims libsvmread.cpp
		mex CFLAGS="\$CFLAGS -std=c99" -largeArrayDims libsvmwrite.cpp
		mex CFLAGS="\$CFLAGS -std=c99" -largeArrayDims svmtrain.cpp ../svm.cpp svm_model_matlab.cpp
		mex CFLAGS="\$CFLAGS -std=c99" -largeArrayDims svmpredict.cpp ../svm.cpp svm_model_matlab.cpp
	end
catch
	fprintf('If make.m fails, please check README about detailed instructions.\n');
end
