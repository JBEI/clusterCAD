import {
  UPDATE_MOLECULE,
  UPDATE_MODULE,
  BEGIN_DOMAIN_SEARCH,
  HANDLE_DOMAIN_SEARCH_RESPONSE,
} from '../actions/actions';

const initialState = {
  sessionID: "",
  jobs: [],
  molecule: "",
  domainSearchObject: {}, // this will need to be replaced with the default configuration
  domainSearchResponse: {},
};

function rootReducer(state = initialState, action) {
  switch (action.type) {
    case UPDATE_MOLECULE:
      return { ...state, molecule: action.payload.molecule };
    case UPDATE_MODULE:
      return { state };
    case BEGIN_DOMAIN_SEARCH:
      return { ...state, domainSearchObject: action.payload.modules };
    case HANDLE_DOMAIN_SEARCH_RESPONSE:
      return { ...state, domainSearchResponse: action.payload.responseObject };
    default:
      return state;
  }
}

export default rootReducer;