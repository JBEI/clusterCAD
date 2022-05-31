// for now this is the only action definition file
// if we get a lot of actions we can break them out 
// by category
// actions must be synchronous, if you need to do something with backend / async put it in middleware folder

// TYPES
export const UPDATE_MOLECULE = 'UPDATE_MOLECULE';
export const BEGIN_DOMAIN_SEARCH = 'BEGIN_DOMAIN_SEARCH';

// CONSTANTS
const NOW = new Date();
  
// ACTION CREATORS

export function updateMolecule(molecule) {
  console.log('action updateMolecule, payload ' + molecule);
  return { 
    type: UPDATE_MOLECULE, 
    payload: {molecule: molecule}
  }
}

export function beginDomainSearch(moduleArray) {
  console.log('action beginDomainSearch, payload ' + moduleArray);
  return {
    type: BEGIN_DOMAIN_SEARCH,
    payload: {
      modules: moduleArray,
      timestamp: NOW,
      // we might need a job id here
    },
  }
}