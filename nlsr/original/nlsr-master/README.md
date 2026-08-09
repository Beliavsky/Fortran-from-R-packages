This repo to allow access to others to working files in nlsr package by J C Nash as of 2022-9-1.

Attempts to alias wrapnlsr() as nlsr() have not been successful. However, explicit copy of code
from wrapnlsr() to nlsr() with name changes gives a sane package. Attempts with 

    nlsr <- wrapnlsr

gave check errors, mainly w.r.t. the .Rd manual files. 

JN 
