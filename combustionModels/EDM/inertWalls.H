
    
    forAll(patches, patchI)
    {
	const fvPatch& p = patches[patchI];

	labelList neighbours(p.size());

	// only do walls
	if (isType<wallFvPatch>(p))
	{
		const polyPatch& pp = p.patch();
		forAll(pp, patchFaceI)
		{
			 scalar index = 0;

			 index = own[pp.start() + patchFaceI];
			 //if (this->mesh().time().value() > 0.0400005) {
			 this->wFuel_[index] = 0.05*this->wFuel_[index];//*0.0;
			 scalar oppPatchFaceI = this->mesh().cells()[index].opposingFaceLabel(pp.start() + patchFaceI,this->mesh().faces());
			 //forAll(this->mesh().cells()[index].faces(),face)
			 //{
				
			 index = own[oppPatchFaceI];
			 //if (this->mesh().time().value() > 0.0400005) {
			 this->wFuel_[index] = 0.001*this->wFuel_[index];//kappa_[index]*0.0;
		}
	}
    }
