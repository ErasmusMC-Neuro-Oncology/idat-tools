#!/bin/bash

idat-tools mix -r 0.5 GSM6379997_203927450093_R01C01_Grn.idat GSM3024450_200392810022_R04C01_Grn.idat /tmp/rmme.idat

idat-tools subtract -r 0.5 /tmp/rmme.idat GSM3024450_200392810022_R04C01_Grn.idat /tmp/rmme_subtracted.idat
