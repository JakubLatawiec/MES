workspace "MES"
    architecture "x64"
    configurations 
    {   
        "Debug",
        "Release" 
    } 

    project "MES"
        location "src"
        kind "ConsoleApp"   
        language "C++"
        cppdialect "C++17"   
        targetdir "bin/%{cfg.buildcfg}" 
        objdir "obj/%{cfg.buildcfg}"
        files 
        { 
            "**.hpp", 
            "**.cpp",
            "**.h"
        } 

        filter "configurations:Debug"
            defines { "DEBUG" }  
            symbols "On" 

        filter "configurations:Release"  
            defines { "NDEBUG" }    
            optimize "On" 