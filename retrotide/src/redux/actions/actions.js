// for now this is the only action definition file
// if we get a lot of actions we can break them out 
// by category
// actions must be synchronous, if you need to do something with backend / async put it in middleware folder

// TYPES
export const UPDATE_MOLECULE = 'UPDATE_MOLECULE';

// CONSTANTS

// ACTION CREATORS

export function updateMolecule(molecule) {
  console.log('action updateMolecule, payload ' + molecule);
  return { type: UPDATE_MOLECULE, 
            payload: {molecule: molecule} }
}

