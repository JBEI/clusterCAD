import {
  UPDATE_MOLECULE,
  BEGIN_DOMAIN_SEARCH,
} from '../actions/actions';

const initialState = {
  sessionID: "",
  jobs: [],
  molecule: "",
};

function rootReducer(state = initialState, action) {
  console.log("REDUCER HIT: Action ->", action);

  switch (action.type) {
    case UPDATE_MOLECULE:
      return { ...state, molecule: action.payload.molecule };
    case BEGIN_DOMAIN_SEARCH:
      return { ...state };
    // case SUBMIT_DESIGN:
    //   return { ...state, }
    default:
      return state;
  }
}

export default rootReducer;