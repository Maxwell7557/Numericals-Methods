function main
    gauss = GaussMethod
    
    %initializationFromUserInput(gauss)
    initializationFromFile(gauss)
    showSystem(gauss)
    applyGaussMethod(gauss)
    showSolutions(gauss)
end