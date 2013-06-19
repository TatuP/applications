maxCo = 
word usedTimeScheme("localEuler");
bool useLocalEuler = true;
if (mesh.schemesDict().readIfPresent("usedTimeScheme", usedTimeScheme))
{
    if (usedTimeScheme == "localEuler")
    {
        Info<< "time scheme: local Euler" << endl;
    }
    else
    {
        FatalErrorIn
        (
            "rhoCentralFoam::readTimeScheme"
        )   << "usedTimeScheme: " << usedTimeScheme
            << " is not a valid choice. "
            << "Use localEuler"
            << abort(FatalError);
    }
} else
{
	Info << "normal transient timescheme" << endl;
}

