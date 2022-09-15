// async actions and calls to the client / backend go here
// redux-thunk is imported toplevel so we can access it here

import {default as client} from 'axios'; // this is the frontend client
// import { connect } from 'react-redux';
import {
  domainSearchResponseHandler, 
  domainSearchResponseErrorHandler,
} from '../actions/actions';

// stub dispatched from the Sequence Search container
// currently only hitting the api for proof of life
const clusterCADSeqSearch = (molecule, token) => {
  console.log('hit api function search with ' + molecule);
  client.get('/api/', {params: {integer: 0}})
    .then((response) => {console.log(response)})
    .catch((error) => {console.log(error.config)});
}

// Domain Search, called from the DomainSearch container
// dispatches with a json object containing all the modules
// gets back the Django response, which we need to dumpinto an iframe or
// onto another page
const clusterCADDomainSearch = (payload, token) => {
  client.post('/api/', 
                {params: {
                  modules: payload,
                }}, 
                {headers: {
                  "X-CSRFTOKEN": token
                }}
              )
  // console.log("async??");
  // dispatch(domainSearchResponseHandler(response));
    .then((response) => {
      console.log("response ***");
      dispatch(domainSearchResponseHandler(response));
    })
    .catch((error) => {
      console.log("error ***");
      dispatch(domainSearchResponseErrorHandler(error));
    }
  );
}

export {clusterCADSeqSearch, clusterCADDomainSearch};
// then import this function in the component, which has no idea it's async
// remember to update reducers if dispatching an action