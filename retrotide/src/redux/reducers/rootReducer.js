import {
  UPDATE_MOLECULE,
  UPDATE_MODULE,
  BEGIN_DOMAIN_SEARCH,
} from '../actions/actions';

const initialState = {
  sessionID: "",
  jobs: [],
  molecule: "",
  domainSearchObject: {}, // this will need to be replaced with the default configuration
};

function rootReducer(state = initialState, action) {
  console.log("REDUCER HIT: Action ->", action);

  switch (action.type) {
    case UPDATE_MOLECULE:
      return { ...state, molecule: action.payload.molecule };
    case UPDATE_MODULE:
      return { state };
    case BEGIN_DOMAIN_SEARCH:
      return { ...state, domainSearchObject: action.payload.modules };
    default:
      return state;
  }
}

export default rootReducer;