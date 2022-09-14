// for now this is the only action definition file
// if we get a lot of actions we can break them out 
// by category
// actions must be synchronous, if you need to do something with backend / async put it in middleware folder

// TYPES
export const UPDATE_MOLECULE = 'UPDATE_MOLECULE';
export const BEGIN_DOMAIN_SEARCH = 'BEGIN_DOMAIN_SEARCH';
export const UPDATE_MODULE = 'UPDATE_MODULE';
export const HANDLE_DOMAIN_SEARCH_RESPONSE = 'HANDLE_DOMAIN_SEARCH_RESPONSE';

// CONSTANTS
const NOW = new Date();
  
// ACTION CREATORS
// not currently in use
export function updateMolecule(molecule) {
  console.log('action updateMolecule, payload ' + molecule);
  return { 
    type: UPDATE_MOLECULE, 
    payload: {molecule: molecule}
  }
}

// used in DomainSearch page to make sure adding domains and changing options is recorded in state
// this way submitted modules are accurate
export function updateModuleDomains(module) {
  console.log('action updateModuleDomains, payload');
  console.log(module);
  return {
    type: UPDATE_MODULE,
    payload: {module: module}
  }
}

// DomainSearch has been dispatched
export function beginDomainSearch(moduleArray) {
  console.log('action beginDomainSearch, payload ');
  console.log(moduleArray);
  return {
    type: BEGIN_DOMAIN_SEARCH,
    payload: {
      modules: moduleArray,
      timestamp: NOW,
      // we might need a job id here
    },
  }
}

// Got a response back from Domain Search
export function domainSearchResponseHandler(response) {
  console.log('action domainSearchResponse');
  console.log(response);
  return {
    type: HANDLE_DOMAIN_SEARCH_RESPONSE,
    payload: {
      reponseObject: response,
      timestamp: NOW,
    },
  }
}

// Got a error back from Domain Search
export function domainSearchResponseErrorHandler(error) {
  console.log('action domainSearchResponseErrorHandler');
  console.log(error);
  return {
    type: HANDLE_DOMAIN_SEARCH_RESPONSE,
    payload: {
      reponseObject: "an error occurred",
      timestamp: NOW,
    },
  }
}