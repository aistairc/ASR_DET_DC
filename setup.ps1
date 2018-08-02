
# Download ksvdbox
$KSVDBOX = "ksvdbox9"
if (-Not (Test-Path "./${KSVDBOX}.zip")) {
    Write-Host "> Download ${KSVDBOX}."
    Invoke-WebRequest -Uri "http://www.cs.technion.ac.il/~ronrubin/Software/${KSVDBOX}.zip" -Outfile ./${KSVDBOX}.zip
}
if (Test-Path "./${KSVDBOX}.zip") {
    Write-Host "> Extract ${KSVDBOX} files."
    Expand-Archive -Path ./${KSVDBOX}.zip -DestinationPath ./${KSVDBOX} -Force
    $copyItems = @{
        Path = "./${KSVDBOX}/private/normrows.m",
               "./${KSVDBOX}/private/normcols.m"
        Destination = "./"
        Force = $true
    }
    Copy-Item @copyItems
    Remove-Item -Path "./${KSVDBOX}" -Recurse -Force
}

# Download CAOLv1.0
if (-Not (Test-Path "./CAOLv1.0.tar.gz")) {
    Write-Host "> Download CAOL."
    Invoke-WebRequest -Uri "http://www.mehrdadya.com/code/CAOLv1.0.tar.gz" -Outfile ./CAOLv1.0.tar.gz
}

# Download CroppedYale
if (-Not (Test-Path "./CroppedYale.zip")) {
    Write-Host "> Download CroppedYale."
    Invoke-WebRequest -Uri "http://vision.ucsd.edu/extyaleb/CroppedYaleBZip/CroppedYale.zip" -Outfile ./CroppedYale.zip
}
if (Test-Path "./CroppedYale.zip") {
    Write-Host "> Expand CroppedYale archive."
    Expand-Archive -Path ./CroppedYale.zip -DestinationPath ./ -Force
    Copy-Item -Path ./CroppedYale/yaleB01/yaleB01_P00A-005E-10.pgm -Destination ./ -Force
}


