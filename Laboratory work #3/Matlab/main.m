function main
    solve = ZeidelAndJakoby;
    
    initializationFromFile(solve);
    showSystem(solve);
    applyJakobyMethod(solve);
    applyZeidelMethod(solve);
end